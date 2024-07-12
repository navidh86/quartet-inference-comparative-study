import os
import subprocess
import pickle, json
import pandas as pd
from IPython.display import display

configurations = ['100gene-100bp', '100gene-1000bp', '1000gene-100bp', '1000gene-1000bp']

REPLICATES = 10

base_folder = os.path.dirname(os.path.dirname((os.path.abspath(__file__))))

def concat_true_gene_trees():
    genes_covered = []
    for config in configurations:
        GENES = int(config.split("-")[0].split("gene")[0])
        if GENES not in genes_covered:
            genes_covered.append(GENES)
        else:
            continue
        for replicate_num in range(1, REPLICATES+1):
            print(f"Replicate {replicate_num}......")
            combined_gene_tree_file = f"{base_folder}/15-taxon/true-genetrees/R{replicate_num}/{GENES}_gt.tre"
            with open(combined_gene_tree_file, "w") as fout:
                for gene in range(1, GENES+1):
                    gene_tree_file = f"{base_folder}/15-taxon/true-genetrees/R{replicate_num}/{gene}/{gene}.tre"

                    with open(gene_tree_file, "r") as fin:
                        tree = fin.readline()
                        # remove all spaces
                        tree = tree.replace(" ", "")

                    fout.write(tree)


def generate_true_gtf_wqrts():
    genes_covered = []
    for config in configurations:
        GENES = int(config.split("-")[0].split("gene")[0])
        if GENES not in genes_covered:
            genes_covered.append(GENES)
        else:
            continue
        for replicate_num in range(1, REPLICATES+1):
            print(f"Replicate {replicate_num}......")
            
            combined_gene_tree_file = f"{base_folder}/15-taxon/true-genetrees/R{replicate_num}/{GENES}_gt.tre"
            output_wqrts_file = f"{base_folder}/15-taxon/Output/true_gt/{GENES}/R{replicate_num}-GTF-true.wqrts"

            cmd = f"python3 generate_wqrts.py {combined_gene_tree_file} {output_wqrts_file}"
            print(cmd)
            os.system(cmd)

def generate_true_dominant_quartets():
    genes_covered = []
    for config in configurations:
        GENES = int(config.split("-")[0].split("gene")[0])
        if GENES not in genes_covered:
            genes_covered.append(GENES)
        else:
            continue
        for replicate_num in range(1, REPLICATES+1):
            input_file = f"{base_folder}/15-taxon/Output/true_gt/{GENES}/R{replicate_num}-GTF-true.wqrts"
            output_file = f"{base_folder}/15-taxon/Output/true_gt/{GENES}/R{replicate_num}-GTF-true.dqrts"
            cmd = f"java -jar {base_folder}/Scripts-2/generateBestWQrts.jar {input_file} {output_file}"
            print(cmd)
            os.system(cmd)

def compute_scores():
    data = {}
    mapping2 = {
    "astral-regular": "ASTRAL",
    'wQFM-GTF': 'wQFM-GTF',
    'wQMC-GTF': 'wQMC-GTF', # (9 July 2024: to compare between wqfm and wqmc on RQ2)
    "wQFM-GTF-bucky-boot": "wQFM-GTF-MB",
    "wQMC-GTF-bucky-boot": "wQMC-GTF-MB",
    'wQFM-SVD-reciprocal': 'wQFM-SVD-Rec',
    'wQMC-SVD-reciprocal': 'wQMC-SVD-Rec',
    'true': 'True tree',
    }
    methods = ['wQFM-GTF', 'astral-regular', 'wQFM-GTF-bucky-boot', 'true', 'wQMC-GTF',
               'wQMC-GTF-bucky-boot', 'wQFM-SVD-reciprocal', 'wQMC-SVD-reciprocal']

    for config in configurations:
        print("Config:", config)
        data[config] = {}
        GENES = int(config.split("-")[0].split("gene")[0])        
        for method in methods:
            mapped_method = mapping2[method]
            data[config][mapped_method] = {}
            total_score_est = 0
            total_score_true = 0
            ratio_score_est = 0
            ratio_score_true = 0
            # true_tree_score = 0
            for replicate_num in range(1, REPLICATES+1):
                est_qrts_file = f"{base_folder}/15-taxon/Output/{config}/R{replicate_num}-GTF.wqrts"
                true_qrts_file = f"{base_folder}/15-taxon/Output/true_gt/{GENES}/R{replicate_num}-GTF-true.wqrts"
                if method != "true":
                    stree_file = f"{base_folder}/Estimated-Species-Trees/15-taxon/{config}/R{replicate_num}-{method}.tre"
                else:
                    stree_file = f"{base_folder}/15-taxon/true-species.tre"

                temp_file = f"{base_folder}/Estimated-Species-Trees/15-taxon/{config}/R{replicate_num}-{method}.temp"
                with open(stree_file, "r") as fin:
                    with open(temp_file, "w") as fout:
                        tree = fin.read()
                        tree = tree.strip()
                        fout.write(tree)


                result_est = subprocess.run(['python3', 'compute_quartet_score.py', est_qrts_file, temp_file, "2"], stdout=subprocess.PIPE).stdout.decode('utf-8')[:-1].split("\t")
                result_true = subprocess.run(['python3', 'compute_quartet_score.py', true_qrts_file, temp_file, "2"], stdout=subprocess.PIPE).stdout.decode('utf-8')[:-1].split("\t")

                print("est", result_est)
                print("true", result_true)
        
                total_score_est += float(result_est[0])
                total_score_true += float(result_true[0])

                ratio_score_est += float(result_est[2])
                ratio_score_true += float(result_true[2])

                os.remove(temp_file)

            print(config, "-->", method)
            print("Estimated gt score:", total_score_est/REPLICATES)
            print("True gt score:", total_score_true/REPLICATES)
            print("\n\n")
            data[config][mapped_method]["est"] = ratio_score_est/REPLICATES
            data[config][mapped_method]["true"] = ratio_score_true/REPLICATES
            data[config][mapped_method]["est_total"] = total_score_est/REPLICATES
            data[config][mapped_method]["true_total"] = total_score_true/REPLICATES


    with open(f"{base_folder}/Results/data/Quartet-Scores/15-taxon/quartet-scores.pkl", "wb") as fp:
        pickle.dump(data, fp)

    with open(f"{base_folder}/Results/data/Quartet-Scores/15-taxon/quartet-scores.json", "w") as fp:
        fp.write(json.dumps(data))
    
    print(data)

def make_df():
    with open(f"{base_folder}/Results/data/Quartet-Scores/15-taxon/quartet-scores.json", "r") as fp:
        results = json.load(fp)

    data = []
    methods = ['wQFM-GTF', 'ASTRAL', 'wQFM-GTF-MB', 'True tree']
    for configuration in configurations:
        for ref in ["est_total", "true_total"]:
            tmp = [configuration, ref]
            for method in methods:
                tmp.append(results[configuration][method][ref])
            data.append(tmp)

    columns = ["Configuration", "Type"] + methods
    df = pd.DataFrame(data, columns=columns)
    display(df)

    result_df_file = f"{base_folder}/Results/data/Quartet-Scores/15-taxon/quartet-scores.xlsx"
    df.to_excel(result_df_file, index=False)

def make_df_2():
    with open(f"{base_folder}/Results/data/Quartet-Scores/15-taxon/quartet-scores.json", "r") as fp:
        results = json.load(fp)

    data = []
    # methods = ['wQFM-GTF', 'ASTRAL', 'wQFM-GTF-MB', 'True tree']
    # methods = ['wQMC-GTF', 'wQFM-GTF', 'ASTRAL', 'wQFM-GTF-MB', 'True tree'] # (9 July 2024: to compare between wqfm and wqmc on RQ2)
    methods = ['wQFM-GTF', 'wQMC-GTF', 'True tree'] # (9 July 2024: to compare between wqfm and wqmc on RQ2)
    for configuration in configurations:
        tmp = [configuration]
        for ref in ["est_total", "true_total"]:
            for method in methods:
                tmp.append(results[configuration][method][ref])
        data.append(tmp)

    columns = ["Configuration"] + [x + "-est" for x in methods] + [x + "-true" for x in methods]
    df = pd.DataFrame(data, columns=columns)
    display(df)

    result_df_file = f"{base_folder}/Results/data/Quartet-Scores/15-taxon/quartet-scores-2.xlsx"
    df.to_excel(result_df_file, index=False)

def make_df_3():
    with open(f"{base_folder}/Results/data/Quartet-Scores/15-taxon/quartet-scores.json", "r") as fp:
        results = json.load(fp)

    data = []
    # methods = ['wQFM-GTF', 'ASTRAL', 'wQFM-GTF-MB', 'True tree']
    # methods = ['wQMC-GTF', 'wQFM-GTF', 'ASTRAL', 'wQFM-GTF-MB', 'True tree'] # (9 July 2024: to compare between wqfm and wqmc on RQ2)
    methods = ['wQFM-GTF', 'wQMC-GTF', 'wQFM-SVD-Rec', 'wQMC-SVD-Rec', 'wQFM-GTF-MB', 'wQMC-GTF-MB', 'True tree'] # (9 July 2024: to compare between wqfm and wqmc on RQ2)
    for configuration in configurations:
        tmp = [configuration]
        # for ref in ["est_total", "true_total"]:
        for ref in ["true_total"]:
            for method in methods:
                tmp.append(results[configuration][method][ref])
        data.append(tmp)

    # columns = ["Configuration"] + [x + "-est" for x in methods] + [x + "-true" for x in methods]
    columns = ["Configuration"] +  [x + "-true" for x in methods]
    df = pd.DataFrame(data, columns=columns)
    display(df)

    result_df_file = f"{base_folder}/Results/data/Quartet-Scores/15-taxon/quartet-scores-3.xlsx"
    df.to_excel(result_df_file, index=False)


# concat_true_gene_trees()
# generate_true_gtf_wqrts()
# generate_true_dominant_quartets() # (9 July 2024: was previously named generate_dominant_quartets())
# compute_scores()
# make_df()
# make_df_2()
make_df_3()