import logging
import pandas as pd
import datetime



def report_extraction(model, df):
    """
    Generate an HTML report summarizing the kcat extraction results.
    
    Parameters:
        model (Model): Metabolic model object
        df (pd.DataFrame): DataFrame containing kcat extraction results

    """
    # Model statistics
    num_model_reactions = len(model.reactions)
    num_model_metabolites = len(model.metabolites)
    num_model_genes = len(model.genes)

    # Kcat extraction statistics
    num_reactions = df['rxn'].nunique()
    num_ec_codes = df['ec_code'].nunique()
    num_kegg_rxn_ids = df['KEGG_rxn_id'].nunique()
    num_ec_kegg_pairs = df[['ec_code', 'KEGG_rxn_id']].drop_duplicates().shape[0]

    # Coverage statistics
    rxn_coverage = 100.0 * num_reactions / num_model_reactions if num_model_reactions else 0

    # Percentage of unique EC codes with at least one KEGG gene
    ec_with_kegg_gene = df.groupby('ec_code')['kegg_genes'].apply(lambda x: any(g for g in x if g)).sum()
    percent_ec_with_kegg_gene = 100.0 * ec_with_kegg_gene / num_ec_codes if num_ec_codes else 0

    html = f"""
    <html>
    <head><title>Kcat Extraction Report</title></head>
    <body>
    <h1>Kcat Extraction Report</h1>
    <h2>Model Overview</h2>
    <ul>
        <li><b>Model name:</b> {model.name}</li>
        <li><b>Number of reactions:</b> {num_model_reactions}</li>
        <li><b>Number of metabolites:</b> {num_model_metabolites}</li>
        <li><b>Number of genes:</b> {num_model_genes}</li>
    </ul>
    <h2>Kcat Extraction Statistics</h2>
    <ul>
        <li><b>Number of reactions with kcat informations:</b> {num_reactions} ({rxn_coverage:.1f}% of model reactions)</li>
        <li><b>Number of unique EC codes:</b> {num_ec_codes}</li>
        <li><b>Number of unique KEGG reaction IDs:</b> {num_kegg_rxn_ids}</li>
        <li><b>Number of unique (EC code, KEGG rxn ID) pairs:</b> {num_ec_kegg_pairs}</li>
        <li><b>Total rows in output:</b> {len(df)}</li>
        <li><b>Percentage of unique EC codes with KEGG genes:</b> {percent_ec_with_kegg_gene:.1f}%</li>
    </ul>
    </body>
    </html>
    """

    # Save report
    report_path = "reports/kcat_summary.html"
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(html)
    logging.info(f"HTML report saved to '{report_path}'")



def report_sabio_rk(df): 
    """
    Generate an HTML report summarizing the kcat matching results from SABIO-RK.

    Parameters:
        df (pd.DataFrame): DataFrame containing kcat matching results. 
    """
    # Gather statistics
    total = len(df)
    matched = df['kcat'].notna().sum()
    match_percent = matched / total * 100 if total > 0 else 0
    score_counts = df['matching_score'].value_counts().sort_index()
    score_percent = (score_counts / total * 100).round(2)

    # Prepare HTML report
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    html = f"""
    <html>
    <head><title>SABIO-RK kcat Matching Report</title></head>
    <body>
        <h1>SABIO-RK kcat Matching Report</h1>
        <p><b>Execution Time:</b> {now}</p>
        <p><b>Total entries:</b> {total}</p>
        <p><b>Matched kcat:</b> {matched} ({match_percent:.2f}%)</p>
        <h2>Matching Score Distribution</h2>
        <table border="1">
            <tr><th>Score</th><th>Count</th><th>Percent</th></tr>
    """
    for score, count in score_counts.items():
        percent = score_percent[score]
        html += f"<tr><td>{score}</td><td>{count}</td><td>{percent:.2f}%</td></tr>"
    html += """
        </table>
    </body>
    </html>
    """

    # Save report
    report_path = "reports/sabio_rk_report.html"
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(html)
    logging.info(f"HTML report saved to '{report_path}'")