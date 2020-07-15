from SPARQLWrapper import SPARQLWrapper, JSON

import pandas as pd    
    
import pickle, hashlib    
    
class QTLSEARCH:
    
    def __init__(self, search, qtls, go_annotations):
        self.qtls = qtls
        self.search = search
        self.go_annotations = go_annotations
        self.p_uniprot_reviewed = 1.00
        self.p_uniprot_unreviewed = 0.95
        self.loss_up_ortholog = 0.85
        self.loss_down_ortholog = 0.85
        self.loss_up_paralog = 0.7225
        self.loss_down_paralog = 0.7225
        #actions        
        print("\033[1m" + "=== GET DATA ===" + "\033[0m")
        self.qtl_gene_roots, self.qtl_gene_protein, self.hog_group_trees, self.hog_group_genes = self.__collect_data()
        print("\033[1m" + "=== COMPUTATIONS ===" + "\033[0m")
        self.__do_computations()
        print("\033[1m" + "=== CREATED QTLSEARCH OBJECT ===" + "\033[0m")
        
    def report(self):                    
        reports = []
        for i in range(0,len(self.qtls)):
            report = []
            for gene in self.qtls[i]:
                if gene in self.qtl_gene_roots.keys():
                    if self.qtl_gene_protein[gene] in self.hog_group_genes[self.qtl_gene_roots[gene]].index :
                        p_initial = self.hog_group_genes[self.qtl_gene_roots[gene]].loc[self.qtl_gene_protein[gene],"p_initial"]
                        p_final = self.hog_group_genes[self.qtl_gene_roots[gene]].loc[self.qtl_gene_protein[gene],"p_final"]
                        report.append([gene, p_initial, p_final, self.qtl_gene_protein[gene]])                                                                    
            df = pd.DataFrame(report)  
            df.columns = ["gene", "p_initial", "p_final", "protein" ]
            df = df.set_index("gene")
            df.sort_values(by=["p_final","p_initial","gene"], ascending=[0, 0, 1], inplace=True)
            reports.append(df)                                          
            return(reports)        
        
    def __collect_data(self):
        #define variables
        hog_group_trees = pd.Series()
        hog_group_genes = pd.Series()
        #root in hog tree for each gene
        qtl_gene_roots = pd.Series()
        qtl_gene_protein = pd.Series()
        gene_p_initial = pd.Series()
        #start collecting
        for qtl in self.qtls:
            for gene in qtl:
                if not(gene in qtl_gene_protein.keys()):
                    p_qtl_initial = 1.0/len(qtl)        
                    #start searching
                    print("\033[1m"+"Search for "+gene+"\033[0m")
                    #first, go up in the hog-tree
                    parent_groups = self.search.get_parent_groups(gene)
                    if len(parent_groups)>0:
                        #define the root of the tree
                        qtl_gene_roots[gene] = parent_groups.index[0]
                        qtl_gene_protein[gene] = parent_groups.loc[parent_groups.index[0]].protein
                        print("- root is "+qtl_gene_roots[gene])
                        if qtl_gene_roots[gene] in hog_group_trees.index:
                            print("- tree already created")
                            if qtl_gene_protein[gene] in gene_p_initial.index:
                                hog_group_genes[qtl_gene_roots[gene]].loc[qtl_gene_protein[gene],"p_initial"] = max(gene_p_initial[qtl_gene_protein[gene]], p_qtl_initial)
                            else:
                                hog_group_genes[qtl_gene_roots[gene]].loc[qtl_gene_protein[gene],"p_initial"] = p_qtl_initial
                        else:    
                            #go down the hog-tree, just to find tree structure
                            hog_group_trees[qtl_gene_roots[gene]] = self.search.get_child_groups(qtl_gene_roots[gene])
                            hog_group_trees[qtl_gene_roots[gene]].loc[:,"p_initial"] = pd.Series(0.0, index=hog_group_trees[qtl_gene_roots[gene]].index)        
                            hog_group_trees[qtl_gene_roots[gene]].loc[:,"p_up"] = pd.Series(0.0, index=hog_group_trees[qtl_gene_roots[gene]].index)        
                            hog_group_trees[qtl_gene_roots[gene]].loc[:,"p_down"] = pd.Series(0.0, index=hog_group_trees[qtl_gene_roots[gene]].index)        
                            print("- tree of groups fetched: "+str(len(hog_group_trees[qtl_gene_roots[gene]])))
                            #go down again, now to find proteins
                            tree_proteins = self.search.get_child_proteins(qtl_gene_roots[gene])
                            tree_proteins_uniprot = self.search.get_child_proteins_uniprot(qtl_gene_roots[gene])
                            print("- proteins within tree fetched: "+str(len(tree_proteins)))
                            print("- uniprot proteins within tree fetched: "+str(len(tree_proteins_uniprot)))
                            #create final list of checked proteins
                            hog_group_genes[qtl_gene_roots[gene]] = tree_proteins
                            hog_group_genes[qtl_gene_roots[gene]].loc[:,"reviewed"] = pd.Series("unknown", index=hog_group_genes[qtl_gene_roots[gene]].index)
                            hog_group_genes[qtl_gene_roots[gene]].loc[:,"p_initial"] = pd.Series(0.0, index=hog_group_genes[qtl_gene_roots[gene]].index)
                            hog_group_genes[qtl_gene_roots[gene]].loc[:,"p_up"] = pd.Series(0.0, index=hog_group_genes[qtl_gene_roots[gene]].index)
                            hog_group_genes[qtl_gene_roots[gene]].loc[:,"p_final"] = pd.Series(0.0, index=hog_group_genes[qtl_gene_roots[gene]].index)
                            #create uniprot list for checking
                            list_proteins = list(tree_proteins_uniprot.index)
                            #divide proteins into chunks (sparql endpoint doesn't allow huge lists)
                            n=50
                            lists_proteins = [list_proteins[i * n:(i + 1) * n] for i in range((len(list_proteins) + n - 1) // n )]          
                            checked_proteins_list = []
                            #check chunks of proteins in uniprot
                            for i in range(0,len(lists_proteins)):
                                checked_proteins_sublist = self.search.check_uniprot_annotations(lists_proteins[i], list(self.go_annotations.index))
                                print("  * check proteins ("+str((i*n)+1)+"-"+str(min((i+1)*n,len(tree_proteins)))+"): "+
                                      str(len(checked_proteins_sublist))+" with required annotation")
                                #collect group labels to give some indication of origin
                                mask = list(tree_proteins_uniprot.loc[checked_proteins_sublist.index].protein.unique())
                                mask = list(set(tree_proteins.index).intersection(mask))
                                sublist_labels = list(tree_proteins.loc[mask].label.dropna().unique())
                                if len(sublist_labels)>0:
                                      print("    -> "+", ".join(sublist_labels))
                                #add checked proteins from chunk to the collective list        
                                checked_proteins_list.append(checked_proteins_sublist)
                            checked_proteins = pd.concat(checked_proteins_list)                 
                            #update reviewed
                            hog_group_genes[qtl_gene_roots[gene]].reviewed = None
                            if qtl_gene_protein[gene] in gene_p_initial.index:
                                hog_group_genes[qtl_gene_roots[gene]].loc[qtl_gene_protein[gene],"p_initial"] = max(gene_p_initial[qtl_gene_protein[gene]], p_qtl_initial)
                            else:
                                hog_group_genes[qtl_gene_roots[gene]].loc[qtl_gene_protein[gene],"p_initial"] = p_qtl_initial
                            if len(checked_proteins)>0:
                                maskUnreviewed = list(set(tree_proteins_uniprot.index).intersection(checked_proteins.loc[checked_proteins.reviewed=="false"].index))
                                maskUnreviewed = list(set(hog_group_genes[qtl_gene_roots[gene]].index).intersection(tree_proteins_uniprot.loc[maskUnreviewed].protein.unique()))
                                maskUnreviewed = hog_group_genes[qtl_gene_roots[gene]].loc[maskUnreviewed]
                                maskReviewed = list(set(tree_proteins_uniprot.index).intersection(checked_proteins.loc[checked_proteins.reviewed=="true"].index))
                                maskReviewed = list(set(hog_group_genes[qtl_gene_roots[gene]].index).intersection(tree_proteins_uniprot.loc[maskReviewed].protein.unique()))
                                maskReviewed = hog_group_genes[qtl_gene_roots[gene]].loc[maskReviewed]
                                hog_group_genes[qtl_gene_roots[gene]].loc[maskUnreviewed.index,"reviewed"] = False
                                hog_group_genes[qtl_gene_roots[gene]].loc[maskReviewed.index,"reviewed"] = True
                                #update probability                        
                                for reviewedProtein in hog_group_genes[qtl_gene_roots[gene]][hog_group_genes[qtl_gene_roots[gene]].reviewed==True].index:
                                    hog_group_genes[qtl_gene_roots[gene]].loc[reviewedProtein,"p_initial"] = max(hog_group_genes[qtl_gene_roots[gene]].loc[reviewedProtein,"p_initial"],self.p_uniprot_reviewed)
                                for unreviewedProtein in hog_group_genes[qtl_gene_roots[gene]][hog_group_genes[qtl_gene_roots[gene]].reviewed==False].index:
                                    hog_group_genes[qtl_gene_roots[gene]].loc[unreviewedProtein,"p_initial"] = max(hog_group_genes[qtl_gene_roots[gene]].loc[unreviewedProtein,"p_initial"],self.p_uniprot_unreviewed)
                            for positiveProbablilityProtein in hog_group_genes[qtl_gene_roots[gene]].loc[hog_group_genes[qtl_gene_roots[gene]]["p_initial"]>0].index:
                                if positiveProbablilityProtein in gene_p_initial.index:
                                    gene_p_initial[positiveProbablilityProtein] = max(hog_group_genes[qtl_gene_roots[gene]].loc[positiveProbablilityProtein,"p_initial"],gene_p_initial[positiveProbablilityProtein])
                                else:        
                                    gene_p_initial[positiveProbablilityProtein] = hog_group_genes[qtl_gene_roots[gene]].loc[positiveProbablilityProtein,"p_initial"]
                            #report results
                            print("- checked "+str(len(list_proteins))+" uniprot proteins: "+str(len(checked_proteins))+" with required annotation")
                            #including details about reviewed status
                            if len(checked_proteins)>0:
                                reviewedList = checked_proteins.groupby("reviewed")
                                if len(reviewedList)>0:
                                    for label,number in reviewedList.size().items():
                                        print("  -> for "+str(number)+" of these "+str(len(list_proteins))+" items, reviewed is "+str(label))                     
                            print("- checked "+str(len(tree_proteins))+" proteins: "+str(len(hog_group_genes[qtl_gene_roots[gene]]["reviewed"].dropna()))+" linked to uniprot with required annotation")
                            if len(hog_group_genes[qtl_gene_roots[gene]])>0:
                                reviewedList = hog_group_genes[qtl_gene_roots[gene]].groupby("reviewed")
                                if len(reviewedList)>0:
                                    for label,number in reviewedList.size().items():
                                        print("  -> for "+str(number)+" of these "+str(len(hog_group_genes[qtl_gene_roots[gene]]))+" items reviewed is marked as "+str(label))                     

                    else:
                        #no hog groups, can't do anything with this...
                        print("- no groups found")

        #update initial probabilities from other genes and qtls    
        for gene in qtl_gene_protein.keys():
            for protein in hog_group_genes[qtl_gene_roots[gene]].keys():
                if protein in gene_p_initial.keys():
                    hog_group_genes[qtl_gene_roots[gene]].loc[protein,"p_initial"] = gene_p_initial[protein];

        return(qtl_gene_roots, qtl_gene_protein, hog_group_trees, hog_group_genes);
    
    
    def __do_computations(self):
        
        def get_p_up(gene, group):
            children = self.hog_group_trees[self.qtl_gene_roots[gene]].loc[self.hog_group_trees[self.qtl_gene_roots[gene]].parent==group]
            proteins = self.hog_group_genes[self.qtl_gene_roots[gene]].loc[self.hog_group_genes[self.qtl_gene_roots[gene]].group==group]
            type = self.hog_group_trees[self.qtl_gene_roots[gene]].loc[group,"type"]
            p = self.hog_group_trees[self.qtl_gene_roots[gene]].loc[group,"p_initial"]
            for protein in proteins.index:
                if type=="ortholog":
                    p +=  self.loss_up_ortholog * self.hog_group_genes[self.qtl_gene_roots[gene]].loc[protein,"p_initial"]
                elif type=="paralog":    
                    p +=  self.loss_up_paralog * self.hog_group_genes[self.qtl_gene_roots[gene]].loc[protein,"p_initial"]
            for child in children.index:
                if type=="ortholog":
                    p +=  self.loss_up_ortholog * get_p_up(gene, child)
                elif type=="paralog":    
                    p +=  self.loss_up_paralog * get_p_up(gene, child)   
            self.hog_group_trees[self.qtl_gene_roots[gene]].loc[group,"p_up"] = p
            #print("Group "+group+" - "+type+": "+str(p))
            return p;

        def set_p_down(gene, group):
            children = self.hog_group_trees[self.qtl_gene_roots[gene]].loc[self.hog_group_trees[self.qtl_gene_roots[gene]].parent==group]
            proteins = self.hog_group_genes[self.qtl_gene_roots[gene]].loc[self.hog_group_genes[self.qtl_gene_roots[gene]].group==group]
            type = self.hog_group_trees[self.qtl_gene_roots[gene]].loc[group,"type"]
            p = self.hog_group_trees[self.qtl_gene_roots[gene]].loc[group,"p_up"]
            parent = self.hog_group_trees[self.qtl_gene_roots[gene]].loc[group,"parent"]
            if not(parent==None):
              parent_type = self.hog_group_trees[self.qtl_gene_roots[gene]].loc[group,"type"]  
              parent_p = self.hog_group_trees[self.qtl_gene_roots[gene]].loc[parent,"p_down"]
              if parent_type=="ortholog":  
                  p = max(p,self.loss_down_ortholog*parent_p)
              elif parent_type=="paralog":
                  p = max(p,self.loss_down_paralog*parent_p)
            self.hog_group_trees[self.qtl_gene_roots[gene]].loc[group,"p_down"] = p
            for protein in proteins.index:  
                p_protein = self.hog_group_genes[self.qtl_gene_roots[gene]].loc[protein,"p_initial"]
                if type=="ortholog":
                    p_protein = max(self.loss_down_ortholog*p,p_protein)
                elif type=="paralog":
                    p_protein = max(self.loss_down_paralog*p,p_protein)    
                self.hog_group_genes[self.qtl_gene_roots[gene]].loc[protein,"p_final"] = p_protein          
            for child in children.index: 
                set_p_down(gene, child)
    
        for qtl in self.qtls:
            for gene in qtl:        
                if gene in self.qtl_gene_roots.keys():
                    print("Compute for "+gene)
                    #reset
                    self.hog_group_trees[self.qtl_gene_roots[gene]].loc[:,"p_up"] = 0.0
                    self.hog_group_trees[self.qtl_gene_roots[gene]].loc[:,"p_down"] = 0.0
                    self.hog_group_genes[self.qtl_gene_roots[gene]].loc[:,"p_final"] = 0.0
                    #go up
                    get_p_up(gene, self.qtl_gene_roots[gene])
                    #go down
                    set_p_down(gene, self.qtl_gene_roots[gene])
                else:
                    print("Skip "+gene)
        

    
class SEARCH:

    def __init__(self, url_pbg, url_oma, url_uniprot):
        #define cache directory
        self.cache = "cache/"
        #define url
        self.url_pbg = url_pbg
        self.url_oma = url_oma
        self.url_uniprot = url_uniprot
        #define sparql
        self.sparql_pbg = SPARQLWrapper(self.url_pbg)
        self.sparql_pbg.setReturnFormat(JSON)
        self.sparql_oma = SPARQLWrapper(self.url_oma)
        self.sparql_oma.setReturnFormat(JSON)
        self.sparql_uniprot = SPARQLWrapper(self.url_uniprot)
        self.sparql_uniprot.setReturnFormat(JSON)
        
    def cache_name(self, method, parameters) :
        key = method+"_"+hashlib.md5(pickle.dumps(parameters)).hexdigest()
        return(key)
    
    
    def get_location(self, id):
        filename = self.cache + self.cache_name("get_location", id)
        try:
            infile = open(filename,"rb")
            new_object = pickle.load(infile)
            infile.close()
            return(new_object)
        except FileNotFoundError:
            file = open("queries/gene_location.sparql", "r") 
            query = file.read()
            file.close()
            self.sparql_pbg.setQuery(query % id)
            # JSON example
            response = self.sparql_pbg.query().convert()
            result = []
            if response["results"]["bindings"]: 
                for item in response["results"]["bindings"]:
                    result.append([
                    item["gene_id"]["value"],
                    item["location"]["value"],
                    item["begin_ref"]["value"],
                    item["begin_pos"]["value"],
                    item["end_ref"]["value"],
                    item["end_pos"]["value"]])
                df = pd.DataFrame(result)  
                df.columns = ["gene_id", "location", "begin_ref", "begin_pos", "end_ref", "end_pos" ]
                df = df.set_index("gene_id")
                df["begin_pos"] = pd.to_numeric(df["begin_pos"])
                df["end_pos"] = pd.to_numeric(df["end_pos"])
                #cache
                outfile = open(filename,"wb")
                pickle.dump(df, outfile)
                outfile.close()
                return df 
            else:
                return pd.DataFrame()   
        
    def compute_interval(self, g1, g2):  
        locations = pd.concat([self.get_location(g1), self.get_location(g2)])
        display(locations[["location"]])
        if(len(locations.index)!=2) :
            print("unexpected number of rows in locations:",len(locations.index))
        elif(locations.iloc[0]['end_pos']>locations.iloc[1]['begin_pos']) & (g1!=g2) :
            print("unexpected order",locations.index[0],"and",locations.index[1])
        else :
            result = []
            if locations.iloc[0]["end_pos"]>locations.iloc[0]["begin_pos"] :
              result.append(["begin", locations.iloc[0]["end_ref"], locations.iloc[0]["end_pos"]])
            else :
              result.append(["begin", locations.iloc[0]["begin_ref"], locations.iloc[0]["begin_pos"]])
            if locations.iloc[1]["begin_pos"]<locations.iloc[1]["end_pos"] :
              result.append(["end", locations.iloc[1]["begin_ref"], locations.iloc[1]["begin_pos"]])
            else :
              result.append(["end", locations.iloc[1]["end_ref"], locations.iloc[1]["end_pos"]])
            df = pd.DataFrame(result)
            df.columns = ["type", "ref", "pos" ]
            df = df.set_index("type")
            return df    
        
    def make_interval(self, ref, start, end):
        result = []
        result.append(["begin",ref,start])
        result.append(["end",ref,end])
        df = pd.DataFrame(result)
        df.columns = ["type", "ref", "pos" ]
        df = df.set_index("type")
        return df    
        
    def interval_genes(self, interval):
        filename = self.cache + self.cache_name("interval_genes", interval)
        try:
            infile = open(filename,"rb")
            new_object = pickle.load(infile)
            infile.close()
            return(new_object)
        except FileNotFoundError:
            file = open("queries/interval_genes.sparql", "r") 
            query = file.read()
            file.close()
            self.sparql_pbg.setQuery(query % {"beginRef" : interval.loc["begin"]["ref"], "beginPos" : interval.loc["begin"]["pos"], "endRef" : interval.loc["end"]["ref"], "endPos" : interval.loc["end"]["pos"]})
            # JSON example
            response = self.sparql_pbg.query().convert()
            result = []
            if response["results"]["bindings"]: 
                for item in response["results"]["bindings"]:
                    row = []
                    row.append(item["gene_id"]["value"])
                    row.append(item["location"]["value"])
                    result.append(row)
                df = pd.DataFrame(result)  
                df.columns = ["gene_id", "location"]
                df = df.set_index("gene_id")
                #cache
                outfile = open(filename,"wb")
                pickle.dump(df, outfile)
                outfile.close()
                return df 
            else:
                return pd.DataFrame()      
        
        
    def get_parent_groups(self, protein):
        filename = self.cache + self.cache_name("get_parent_groups", protein)
        try:
            infile = open(filename,"rb")
            new_object = pickle.load(infile)
            infile.close()
            return(new_object)
        except FileNotFoundError:
            file = open("queries/parent_groups.sparql", "r") 
            query = file.read()
            file.close()
            self.sparql_oma.setQuery(query % {"protein" : "\""+protein+"\""})
            # JSON example
            response = self.sparql_oma.query().convert()
            result = []
            if response["results"]["bindings"]: 
                for item in response["results"]["bindings"]:
                    result.append([
                    item["level"]["value"],
                    item["group"]["value"],
                    item["type"]["value"],
                    item["protein"]["value"]])
                df = pd.DataFrame(result)  
                df.columns = ["level", "group", "type", "protein" ]
                df.drop_duplicates(subset ="group", keep = "first", inplace = True) 
                df = df.set_index("group")
                df["level"] = pd.to_numeric(df["level"])
                #cache
                outfile = open(filename,"wb")
                pickle.dump(df, outfile)
                outfile.close()
                return df 
            else:
                return pd.DataFrame()  
        
    #get only group paths that end with a uniprot annotated protein
    def get_child_groups(self, parent):
        filename = self.cache + self.cache_name("get_child_groups", parent)
        try:
            infile = open(filename,"rb")
            new_object = pickle.load(infile)
            infile.close()
            return(new_object)
        except FileNotFoundError:
            file = open("queries/child_groups.sparql", "r") 
            query = file.read()
            file.close()
            self.sparql_oma.setQuery(query % {"parent" : "<"+parent+">"})
            # JSON example
            response = self.sparql_oma.query().convert()
            result = []
            if response["results"]["bindings"]: 
                for item in response["results"]["bindings"]:
                    row = [
                      item["group"]["value"],
                      item["type"]["value"]  
                    ]
                    if "parent" in item.keys() :
                        row.append(item["parent"]["value"])
                    else:
                        row.append(None)
                    if "parent_type" in item.keys() :
                        row.append(item["parent_type"]["value"])
                    else:
                        row.append(None) 
                    if "label" in item.keys() :
                        row.append(item["label"]["value"])
                    else:
                        row.append(None)
                    if "parent_label" in item.keys():
                        row.append(item["parent_label"]["value"])
                    else:
                        row.append(None)    
                    result.append(row)                
                df = pd.DataFrame(result)  
                df.columns = ["group", "type", "parent", "parent_type", "label", "parent_label" ]
                df.drop_duplicates(subset ="group", keep = "first", inplace = True) 
                df = df.set_index("group")
                #cache
                outfile = open(filename,"wb")
                pickle.dump(df, outfile)
                outfile.close()
                return df 
            else:
                return pd.DataFrame()    
        
    #get uniprot annotated proteins and their group
    def get_child_proteins_uniprot(self, parent):
        filename = self.cache + self.cache_name("get_child_proteins_uniprot", parent)
        try:
            infile = open(filename,"rb")
            new_object = pickle.load(infile)
            infile.close()
            return(new_object)
        except FileNotFoundError:
            file = open("queries/child_proteins_uniprot.sparql", "r") 
            query = file.read()
            file.close()
            self.sparql_oma.setQuery(query % {"parent" : "<"+parent+">"})
            # JSON example
            response = self.sparql_oma.query().convert()
            result = []
            if response["results"]["bindings"]: 
                for item in response["results"]["bindings"]:
                    row = [
                      item["group"]["value"],                  
                      item["protein"]["value"]
                    ] 
                    if "protein_uniprot" in item.keys() :
                        row.append(item["protein_uniprot"]["value"])
                    else:
                        row.append(None)
                    if "group_label" in item.keys() :
                        row.append(item["group_label"]["value"])
                    else:
                        row.append(None)
                    result.append(row)                
                df = pd.DataFrame(result)  
                df.columns = ["group", "protein", "uniprot", "label" ]
                df.drop_duplicates(subset ="uniprot", keep = "first", inplace = True) 
                df = df.set_index("uniprot")
                #cache
                outfile = open(filename,"wb")
                pickle.dump(df, outfile)
                outfile.close()
                return df 
            else:
                return pd.DataFrame()    
        
    #get proteins and their group
    def get_child_proteins(self, parent):
        filename = self.cache + self.cache_name("get_child_proteins", parent)
        try:
            infile = open(filename,"rb")
            new_object = pickle.load(infile)
            infile.close()
            return(new_object)
        except FileNotFoundError:
            file = open("queries/child_proteins.sparql", "r") 
            query = file.read()
            file.close()
            self.sparql_oma.setQuery(query % {"parent" : "<"+parent+">"})
            # JSON example
            response = self.sparql_oma.query().convert()
            result = []
            if response["results"]["bindings"]: 
                for item in response["results"]["bindings"]:
                    row = [
                      item["group"]["value"],                  
                      item["protein"]["value"]
                    ] 
                    if "group_label" in item.keys() :
                        row.append(item["group_label"]["value"])
                    else:
                        row.append(None)
                    result.append(row)                
                df = pd.DataFrame(result)  
                df.columns = ["group", "protein", "label" ]
                df.drop_duplicates(subset ="protein", keep = "first", inplace = True) 
                df = df.set_index("protein")
                #cache
                outfile = open(filename,"wb")
                pickle.dump(df, outfile)
                outfile.close()
                return df 
            else:
                return pd.DataFrame()    
        
        
    #get child annotations
    def get_child_annotations(self, go_annotation):
        filename = self.cache + self.cache_name("get_child_annotations", go_annotation)
        try:
            infile = open(filename,"rb")
            new_object = pickle.load(infile)
            infile.close()
            return(new_object)
        except FileNotFoundError:
            file = open("queries/child_annotations.sparql", "r") 
            query = file.read()
            file.close()
            self.sparql_pbg.setQuery(query % {"annotation" : go_annotation})
            # JSON example
            response = self.sparql_pbg.query().convert()
            result = []
            if response["results"]["bindings"]: 
                for item in response["results"]["bindings"]:
                    row = [
                      item["go_annotation"]["value"],
                      item["label"]["value"]
                    ]                
                    result.append(row)                
                df = pd.DataFrame(result)  
                df.columns = ["go_annotation", "label" ]
                df.drop_duplicates(subset ="go_annotation", keep = "first", inplace = True) 
                df = df.set_index("go_annotation")
                #cache
                outfile = open(filename,"wb")
                pickle.dump(df, outfile)
                outfile.close()
                return df 
            else:
                return pd.DataFrame()         
        
    #check uniprot for proteins and GO annotations
    def check_uniprot_annotations(self, uniprot_proteins, go_annotations):
        filename = self.cache + self.cache_name("get_uniprot_annotations", [uniprot_proteins, go_annotations])
        try:
            infile = open(filename,"rb")
            new_object = pickle.load(infile)
            infile.close()
            return(new_object)
        except FileNotFoundError:
            file = open("queries/check_uniprot_annotations.sparql", "r") 
            query = file.read()
            file.close()
            self.sparql_uniprot.setQuery(query % {"proteins" : "<"+">,<".join(uniprot_proteins)+">", "annotations": "<"+">,<".join(go_annotations)+">"})
            # JSON example
            response = self.sparql_uniprot.query().convert()
            result = []
            if response["results"]["bindings"]: 
                for item in response["results"]["bindings"]:
                    row = [
                      item["uniprot"]["value"],
                      item["reviewed"]["value"]
                    ]                
                    result.append(row)                
                df = pd.DataFrame(result)  
                df.columns = ["uniprot", "reviewed" ]
                #df["reviewed"] = df["reviewed"].astype("bool")
                df.drop_duplicates(subset ="uniprot", keep = "first", inplace = True) 
                df = df.set_index("uniprot")
                #cache
                outfile = open(filename,"wb")
                pickle.dump(df, outfile)
                outfile.close()
                return df 
            else:
                return pd.DataFrame()
        
    
    

    