import os
import io
import base64
import logging
import datetime
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogFormatter
from io import BytesIO


def report_extraction(model, df, report_statistics):
    """
    Generate a detailed HTML report summarizing the kcat extraction results
    with separate tables for extracted data and EC issues, and visualizations.
    """
    # Model statistics
    nb_model_reactions = len(model.reactions)
    nb_model_metabolites = len(model.metabolites)
    nb_model_genes = len(model.genes)
    unique_ec_codes = []
    for rxn in model.reactions:
        ec_code = rxn.annotation.get('ec-code')
        if ec_code:
            if isinstance(ec_code, str):
                ec_code = [ec_code.strip()]
            elif isinstance(ec_code, list):
                ec_code = [x.strip() for x in ec_code if x.strip()]
            else:
                ec_code = []
            unique_ec_codes.extend(ec_code)
    nb_model_ec_codes = len(set(unique_ec_codes))

    # Kcat statistics
    nb_reactions = df['rxn'].nunique()
    nb_ec_codes = df['ec_code'].nunique()

    nb_ec_codes_transferred = report_statistics.get('transferred_ec_codes', 0)
    nb_ec_codes_incomplete = report_statistics.get('incomplete_ec_codes', 0)
    # nb_ec_codes_no_genes = report_statistics.get('no_genes_ec_codes', 0)
    nb_reactions_dropped = report_statistics.get('nb_of_reactions_due_to_unconsistent_ec', 0)
    nb_lines_dropped = report_statistics.get('nb_of_lines_dropped_due_to_unconsistent_ec', 0)


    rxn_coverage = 100.0 * nb_reactions / nb_model_reactions if nb_model_reactions else 0
    percent_ec_retrieved = 100.0 * nb_ec_codes / nb_model_ec_codes if nb_model_ec_codes else 0

    # Pie Chart
    pie_chart = {
        "Retrieved": nb_ec_codes,
        "Transferred": nb_ec_codes_transferred,
        "Incomplete": nb_ec_codes_incomplete
    }
    
    pie_chart = {k: v for k, v in pie_chart.items() if v > 0} # Remove zero values

    fig, ax = plt.subplots(figsize=(7, 7))
    wedges, texts, autotexts = ax.pie(
        pie_chart.values(),
        labels=None,
        autopct='%1.1f%%',
        startangle=90,
        colors=["#2ecc71", "#e67e22", "#c0392b"],
        textprops={'fontsize': 16}
    )
    ax.axis('equal')
    ax.legend(
        wedges,
        pie_chart.keys(),
        title="EC",
        loc='lower center',
        bbox_to_anchor=(0.5, -0.2),
        ncol=3,
        frameon=False,
        fontsize=16,           
        title_fontsize=18      
    )

    for text in texts:
        text.set_fontsize(16)
    for autotext in autotexts:
        autotext.set_fontsize(16)
    pie_buffer = io.BytesIO()
    plt.savefig(pie_buffer, format='png', bbox_inches='tight')
    plt.close(fig)
    pie_base64 = base64.b64encode(pie_buffer.getvalue()).decode('utf-8')

    # Time
    generated_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Html report
    html = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Extract kcat Report</title>
        <style>
            body {{
                font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
                background-color: #f4f6f9;
                margin: 0;
                padding: 0;
                color: #333;
            }}
            header {{
                background: linear-gradient(90deg, #2c3e50, #2980b9);
                color: #fff;
                padding: 20px;
                text-align: center;
                box-shadow: 0 2px 6px rgba(0,0,0,0.1);
            }}
            header h1 {{
                margin: 0;
                font-size: 2rem;
            }}
            .container {{
                max-width: 1100px;
                margin: 30px auto;
                padding: 20px;
            }}
            .card {{
                background: #fff;
                border-radius: 12px;
                padding: 20px;
                margin-bottom: 20px;
                box-shadow: 0 2px 8px rgba(0,0,0,0.05);
            }}
            .card h2 {{
                margin-top: 0;
                color: #2980b9;
                border-bottom: 2px solid #e6e6e6;
                padding-bottom: 10px;
                font-size: 1.5rem;
            }}
            .stats-grid {{
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
                gap: 15px;
                margin-top: 15px;
            }}
            .stat-box {{
                background: #f9fafc;
                border-radius: 8px;
                padding: 15px;
                text-align: center;
                border: 1px solid #e2e2e2;
            }}
            table {{
                width: 100%;
                border-collapse: collapse;
                margin-top: 20px;
                font-size: 0.95rem;
            }}
            table th, table td {{
                border: 1px solid #ddd;
                padding: 10px;
                text-align: left;
            }}
            table th {{
                background-color: #2980b9;
                color: #fff;
            }}
            table tr:nth-child(even) {{
                background-color: #f2f2f2;
            }}
            .progress {{
                background-color: #ddd;
                border-radius: 10px;
                overflow: hidden;
                height: 18px;
                width: 100%;
                margin-top: 5px;
            }}
            .progress-bar {{
                background-color: #27ae60;
                height: 100%;
                text-align: right;
                padding-right: 5px;
                color: white;
                font-size: 0.8rem;
                line-height: 18px;
            }}
            footer {{
                text-align: center;
                font-size: 0.9rem;
                color: #777;
                padding: 15px;
                margin-top: 20px;
                border-top: 1px solid #ddd;
            }}
        </style>
    </head>
    <body>
        <header>
            <h1>Extract k<sub>cat</sub> Report</h1>
            <p>Generated on {generated_time}</p>
        </header>

        <div class="container">
            <!-- Model Overview -->
            <div class="card">
                <h2>Model Overview</h2>
                <div class="stats-grid">
                    <div class="stat-box">
                        <h3>{model.id}</h3>
                        <p>Model ID</p>
                    </div>
                    <div class="stat-box">
                        <h3>{nb_model_reactions}</h3>
                        <p>Reactions</p>
                    </div>
                    <div class="stat-box">
                        <h3>{nb_model_metabolites}</h3>
                        <p>Metabolites</p>
                    </div>
                    <div class="stat-box">
                        <h3>{nb_model_genes}</h3>
                        <p>Genes</p>
                    </div>
                </div>
            </div>

            <!-- kcat Extraction Table -->
            <div class="card">
                <h2>k<sub>cat</sub> Extraction Statistics</h2>
                <table>
                    <tr>
                        <th>Metric</th>
                        <th>Value</th>
                        <th>Visualization</th>
                    </tr>
                    <tr>
                        <td>Reactions with EC info</td>
                        <td>{nb_reactions} ({rxn_coverage:.1f}%)</td>
                        <td>
                            <div class="progress">
                                <div class="progress-bar" style="width:{rxn_coverage}%;"></div>
                            </div>
                        </td>
                    </tr>
                    <tr>
                        <td>EC codes found in KEGG</td>
                        <td>{nb_ec_codes} ({percent_ec_retrieved:.1f}%)</td>
                        <td>
                            <div class="progress">
                                <div class="progress-bar" style="width:{percent_ec_retrieved}%;"></div>
                            </div>
                        </td>
                    </tr>
                    <tr>
                        <td>Total rows in output</td>
                        <td>{len(df)}</td>
                        <td>-</td>
                    </tr>
                </table>
            </div>

            <!-- EC Issues Table -->
            <div class="card">
                <h2>Issues in EC Assignment</h2>
                <table>
                    <tr>
                        <th>Cases</th>
                        <th>Count</th>
                    </tr>
                    <tr>
                        <td>Transferred EC codes</td>
                        <td>{nb_ec_codes_transferred}</td>
                    </tr>
                    <tr>
                        <td>Incomplete EC codes</td>
                        <td>{nb_ec_codes_incomplete}</td>
                    </tr>
                    <tr>
                        <td>Number of reactions dropped due to inconsistent EC codes</td>
                        <td>{nb_reactions_dropped}</td>
                    </tr>
                    <tr>
                        <td>Number of k<sub>cat</sub> values dropped due to inconsistent EC codes</td>
                        <td>{nb_lines_dropped}</td>
                    </tr>
                </table>
            </div>

            <!-- Pie Chart Section -->
            <div class="card">
                <h2>EC Distribution</h2>
                <img src="data:image/png;base64,{pie_base64}" alt="EC Pie Chart" style="display:block;margin:20px auto;max-width:350px;">
            </div>
        </div>

        <footer>WILDkCAT</footer>
    </body>
    </html>
    """

    # Save report
    os.makedirs("reports", exist_ok=True)
    report_path = "reports/extract_kcat_report.html"
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(html)
    logging.info(f"HTML report saved to '{report_path}'")


def report_retrieval(df):
    """
    Generate a styled HTML report summarizing the kcat matching results,
    including kcat value distribution and matching score repartition.

    Parameters:
        df (pd.DataFrame): Must contain ['kcat', 'matching_score', ...].
    """
    # Ensure numeric kcat values to avoid TypeError on comparisons
    kcat_values = pd.to_numeric(df['kcat'], errors='coerce').dropna()

    # Only use scores present in the data
    present_scores = sorted(df['matching_score'].dropna().unique())
    score_counts = df['matching_score'].value_counts().reindex(present_scores, fill_value=0)
    total = len(df)
    matched = len(kcat_values)
    match_percent = matched / total * 100 if total else 0
    score_percent = (score_counts / total * 100).round(2) if total else pd.Series(0, index=present_scores)

    # Distinct colors for each score (up to 12, then cycle)
    # Gradient colors from green (best score) to red (worst score)
    distinct_colors = [
        "#27ae60",
        "#43b76e",
        "#60c07c",
        "#7cc98a",
        "#98d298",
        "#b5dbb6",
        "#d1e4c4",
        "#f1e9b6",
        "#f7d97c",
        "#f9c74f",
        "#f8961e",
        "#f3722c",
        "#e67e22",
        "#e74c3c",
        "#c0392b",
        "#a93226",
        "#922b21",
        "#7b241c"
    ]

    def score_color(score):
        idx = present_scores.index(score)
        return distinct_colors[idx % len(distinct_colors)]

    generated_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Histogram with stacked bars for scores
    kcat_hist_base64 = ""
    if not kcat_values.empty:
        min_exp = int(np.floor(np.log10(max(1e-6, kcat_values.min()))))
        max_exp = int(np.ceil(np.log10(kcat_values.max())))
        bins = np.logspace(min_exp, max_exp, num=40)
        fig, ax = plt.subplots(figsize=(10, 6))
        # Stacked histogram by score
        hist_data = [pd.to_numeric(df[df['matching_score'] == score]['kcat'], errors='coerce').dropna() for score in present_scores]
        ax.hist(hist_data, bins=bins, stacked=True, color=[score_color(s) for s in present_scores], label=[f"Score {s}" for s in present_scores], edgecolor='white')
        ax.set_xscale('log')
        ax.set_xlabel("kcat", fontsize=14)
        ax.set_ylabel("Count", fontsize=14)
        ax.tick_params(axis='both', labelsize=12)
        ax.set_xlim([10**min_exp / 1.5, 10**max_exp * 1.5])
        ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
        ax.legend(title="Matching Score", fontsize=12)
        plt.tight_layout()
        buf = io.BytesIO()
        plt.savefig(buf, format='png', bbox_inches='tight')
        plt.close(fig)
        kcat_hist_base64 = base64.b64encode(buf.getvalue()).decode('utf-8')

    # HTML start
    html = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Retrieve kcat Report</title>
        <style>
            body {{
                font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
                background-color: #f4f6f9;
                margin: 0;
                padding: 0;
                color: #333;
            }}
            header {{
                background: linear-gradient(90deg, #2c3e50, #2980b9);
                color: #fff;
                padding: 20px;
                text-align: center;
                box-shadow: 0 2px 6px rgba(0,0,0,0.1);
            }}
            header h1 {{
                margin: 0;
                font-size: 2rem;
            }}
            .container {{
                max-width: 1000px;
                margin: 30px auto;
                padding: 20px;
            }}
            .card {{
                background: #fff;
                border-radius: 12px;
                padding: 20px;
                margin-bottom: 20px;
                box-shadow: 0 2px 8px rgba(0,0,0,0.05);
            }}
            .card h2 {{
                margin-top: 0;
                color: #2980b9;
                border-bottom: 2px solid #e6e6e6;
                padding-bottom: 10px;
                font-size: 1.5rem;
            }}
            .stats-grid {{
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
                gap: 15px;
                margin-top: 15px;
            }}
            .stat-box {{
                background: #f9fafc;
                border-radius: 8px;
                padding: 15px;
                text-align: center;
                border: 1px solid #e2e2e2;
            }}
            .stat-box h3 {{
                margin: 0;
                font-size: 1.3rem;
                color: #2c3e50;
            }}
            .stat-box p {{
                margin: 5px 0 0;
                color: #666;
            }}
            table {{
                width: 100%;
                border-collapse: collapse;
                margin-top: 20px;
                font-size: 0.95rem;
            }}
            table th, table td {{
                border: 1px solid #ddd;
                padding: 10px;
                text-align: left;
            }}
            table th {{
                background-color: #2980b9;
                color: #fff;
            }}
            table tr:nth-child(even) {{
                background-color: #f2f2f2;
            }}
            .progress-stacked {{
                display: flex;
                height: 18px;
                border-radius: 10px;
                overflow: hidden;
                background-color: #ddd;
                font-size: 0.75rem;
                line-height: 18px;
                color: white;
                text-shadow: 0 1px 1px rgba(0,0,0,0.2);
                margin-bottom: 10px;
            }}
            .progress-bar {{
                display: flex;
                align-items: center;
                justify-content: center;
                height: 100%;
                white-space: nowrap;
                overflow: hidden;
            }}
            .legend {{
                display: flex;
                flex-wrap: wrap;
                gap: 10px;
                font-size: 0.85rem;
                margin-top: 5px;
            }}
            .legend-item {{
                display: flex;
                align-items: center;
                gap: 5px;
            }}
            .legend-color {{
                width: 14px;
                height: 14px;
                border-radius: 3px;
                border: 1px solid #aaa;
            }}
            .img-section {{
                display: flex;
                flex-wrap: wrap;
                gap: 30px;
                justify-content: center;
                align-items: flex-start;
                margin-top: 20px;
            }}
            footer {{
                text-align: center;
                font-size: 0.9rem;
                color: #777;
                padding: 15px;
                margin-top: 20px;
                border-top: 1px solid #ddd;
            }}
        </style>
    </head>
    <body>
        <header>
            <h1>Retrieve k<sub>cat</sub> Report</h1>
            <p>Generated on {generated_time}</p>
        </header>

        <div class="container">
            <div class="card">
                <h2>Overview</h2>
                <div class="stats-grid">
                    <div class="stat-box">
                        <h3>{total}</h3>
                        <p>Total Entries</p>
                    </div>
                    <div class="stat-box">
                        <h3>{matched}</h3>
                        <p>Matched k<sub>cat</sub> ({match_percent:.2f}%)</p>
                    </div>
                </div>
            </div>

            <div class="card">
                <h2>Matching Score Distribution</h2>
                <div class="progress-stacked">
    """

    # Add progress bars only for present scores
    for score in present_scores:
        percent = score_percent.get(score, 0)
        if percent > 0:
            html += f'<div class="progress-bar" style="width:{percent}%;background:{score_color(score)};" title="Score {score}: {percent:.2f}%"></div>'

    html += """
            </div>
            <div class="legend">
    """

    # Add legend only for present scores
    for score in present_scores:
        html += f'<div class="legend-item"><div class="legend-color" style="background:{score_color(score)};"></div> Score {score}</div>'

    html += """
            </div>
            <table>
                <tr>
                    <th>Score</th>
                    <th>Count</th>
                    <th>Percent</th>
                </tr>
    """

    # Table rows only for present scores
    for score in present_scores:
        html += f'<tr><td>{score}</td><td>{score_counts[score]}</td><td>{score_percent[score]:.2f}%</td></tr>'

    html += """
            </table>
        </div>
    """

    # Histogram section (stacked by score)
    html += """
        <div class="card">
            <h2>Distribution of k<sub>cat</sub> values (Stacked by Matching Score)</h2>
            <div class="img-section">
    """
    if kcat_hist_base64:
        html += f'<img src="data:image/png;base64,{kcat_hist_base64}" alt="k<sub>cat</sub> Distribution">'
    html += """
            </div>
        </div>
    """

    # Metadata section
    html += """
            <div class="card">
                <h2>Matching Score</h2>
                <p>
                    The matching score evaluates how well a candidate k<sub>cat</sub> entry fits the query enzyme and conditions. 
                    A lower score indicates a better match (0 = best possible, 15 = no match).
                </p>
                <h3>Scoring process:</h3>
                <ul>
                    <li><b>Catalytic enzyme:</b> Check if the reported enzyme matches the expected catalytic enzyme(s).</li>
                    <li><b>Organism:</b> Penalize mismatches between the source organism and the target organism.</li>
                    <li><b>Enzyme variant:</b> Exclude or penalize mutant/engineered variants (wildtype preferred).</li>
                    <li><b>pH:</b> Check whether the reported pH is consistent with the desired experimental range.</li>
                    <li><b>Substrate:</b> Verify substrate compatibility with the catalytic reaction.</li>
                    <li><b>Temperature:</b> Penalize deviations from the target temperature; 
                        if possible, adjust kcat values using the Arrhenius equation.</li>
                </ul>

                <h3>Score breakdown (default penalties):</h3>
                <table border="1" cellpadding="6" cellspacing="0" style="border-collapse: collapse; text-align: left;">
                    <tr>
                        <th>Criterion</th>
                        <th>Penalty</th>
                    </tr>
                    <tr>
                        <td>Substrate mismatch</td>
                        <td>+3</td>
                    </tr>
                    <tr>
                        <td>Catalytic enzyme mismatch</td>
                        <td>+2</td>
                    </tr>
                    <tr>
                        <td>Organism mismatch</td>
                        <td>+2</td>
                    </tr>
                    <tr>
                        <td>pH unknown</td>
                        <td>+1</td>
                    </tr>
                    <tr>
                        <td>pH out of range</td>
                        <td>+2</td>
                    </tr>
                    <tr>
                        <td>Temperature unknown</td>
                        <td>+1</td>
                    </tr>
                    <tr>
                        <td>Temperature out of range</td>
                        <td>+2</td>
                    </tr>
                    <tr>
                        <td>Enzyme variant unknown</td>
                        <td>+1</td>
                    </tr>
                </table>

                <p>
                    Candidates are then ranked by:
                    <ol>
                        <li>Lowest total score</li>
                        <li>Highest sequence identity percentage to the target enzyme</li>
                        <li>Adjusted k<sub>cat</sub> value (favoring the smallest value by default)</li>
                    </ol>
                </p>
                <p>
                    The best candidate is the one with the lowest score after these checks. 
                    If multiple candidates tie on score, sequence identity and k<sub>cat</sub> values break the tie.
                </p>
            </div>
        </div>

        <footer>WILDkCAT</footer>
    </body>
    </html>
    """

    # Save HTML
    os.makedirs("reports", exist_ok=True)
    report_path = f"reports/retrieve_kcat_report.html"
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(html)

    logging.info(f"HTML report saved to '{report_path}'")


def report_prediction_input(catapro_df, report_statistics): 
    # CataPro Statistics 
    total_catapro_entries = len(catapro_df)

    # Report Statistics
    rxn_covered = report_statistics['rxn_covered']
    cofactors_covered = report_statistics['cofactor_identified']
    no_catalytic = report_statistics['no_catalytic']
    kegg_missing = report_statistics['kegg_no_matching']

    total_rxn = rxn_covered + kegg_missing + no_catalytic
    rxn_coverage = (rxn_covered / total_rxn * 100) if total_rxn > 0 else 0

    # Time
    generated_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Html report
    html = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Predict kcat Report - Part 1</title>
        <style>
            body {{
                font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
                background-color: #f4f6f9;
                margin: 0;
                padding: 0;
                color: #333;
            }}
            header {{
                background: linear-gradient(90deg, #2c3e50, #2980b9);
                color: #fff;
                padding: 20px;
                text-align: center;
                box-shadow: 0 2px 6px rgba(0,0,0,0.1);
            }}
            header h1 {{
                margin: 0;
                font-size: 2rem;
            }}
            .container {{
                max-width: 1100px;
                margin: 30px auto;
                padding: 20px;
            }}
            .card {{
                background: #fff;
                border-radius: 12px;
                padding: 20px;
                margin-bottom: 20px;
                box-shadow: 0 2px 8px rgba(0,0,0,0.05);
            }}
            .card h2 {{
                margin-top: 0;
                color: #2980b9;
                border-bottom: 2px solid #e6e6e6;
                padding-bottom: 10px;
                font-size: 1.5rem;
            }}
            .stats-grid {{
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
                gap: 15px;
                margin-top: 15px;
            }}
            .stat-box {{
                background: #f9fafc;
                border-radius: 8px;
                padding: 15px;
                text-align: center;
                border: 1px solid #e2e2e2;
            }}
            table {{
                width: 100%;
                border-collapse: collapse;
                margin-top: 20px;
                font-size: 0.95rem;
            }}
            table th, table td {{
                border: 1px solid #ddd;
                padding: 10px;
                text-align: left;
            }}
            table th {{
                background-color: #2980b9;
                color: #fff;
            }}
            table tr:nth-child(even) {{
                background-color: #f2f2f2;
            }}
            .progress {{
                background-color: #ddd;
                border-radius: 10px;
                overflow: hidden;
                height: 18px;
                width: 100%;
                margin-top: 5px;
            }}
            .progress-bar {{
                background-color: #27ae60;
                height: 100%;
                text-align: right;
                padding-right: 5px;
                color: white;
                font-size: 0.8rem;
                line-height: 18px;
            }}
            footer {{
                text-align: center;
                font-size: 0.9rem;
                color: #777;
                padding: 15px;
                margin-top: 20px;
                border-top: 1px solid #ddd;
            }}
        </style>
    </head>
    <body>
        <header>
            <h1>Predict k<sub>cat</sub> Report - Part 1</h1>
            <p>Generated on {generated_time}</p>
        </header>

        <div class="container">
            <!-- CataPro Overview -->
            <div class="card">
                <h2>Overview</h2>
                <div class="stats-grid">
                    <div class="stat-box">
                        <h3>{total_rxn}</h3>
                        <p>Total k<sub>cat</sub> values</p>
                    </div>
                    <div class="stat-box">
                        <h3>{rxn_covered}</h3>
                        <p>k<sub>cat</sub> to be predicted ({rxn_coverage:.2f}%)</p>
                    </div>
                </div>
            </div>

            <!-- Prediction kcat Table -->
            <div class="card">
                <h2>k<sub>cat</sub> Prediction Statistics</h2>
                <table>
                    <tr>
                        <th>Metric</th>
                        <th>Value</th>
                    </tr>
                    <tr>
                        <td>Total of entries in CataPro input file</td>
                        <td>{total_catapro_entries}</td>
                    </tr>
                    <tr>
                        <td>Number of cofactor identified</td>
                        <td>{cofactors_covered}</td>
                    </tr>
                    <tr>
                        <td>Entries with no catalytic activity identified</td>
                        <td>{no_catalytic}</td>
                    </tr>
                    <tr>
                        <td>Entries with missing KEGG IDs</td>
                        <td>{kegg_missing}</td>
                    </tr>
                </table>
            </div>
            
        <footer>WILDkCAT</footer>
    </body>
    </html>
    """

    # Save report
    os.makedirs("reports", exist_ok=True)
    report_path = "reports/prediction_kcat_report.html"
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(html)
    logging.info(f"HTML report saved to '{report_path}'")


def report_final(final_df):
    """
    Generate a full HTML report summarizing predicted vs. source kcat values.
    """

    df = final_df.copy()
    df["kcat_db"] = df["kcat_db"].fillna("Unknown")
    generated_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Utility to convert matplotlib figures to base64 <img>
    def fig_to_base64(fig):
        buf = BytesIO()
        fig.savefig(buf, format="png", bbox_inches="tight")
        buf.seek(0)
        encoded = base64.b64encode(buf.read()).decode("utf-8")
        plt.close(fig)
        return f'<div class="plot-container"><img src="data:image/png;base64,{encoded}"></div>'

    # === 1. Distribution plots ===
    def plot_kcat_distribution(column_name, title):
        df[column_name] = pd.to_numeric(df[column_name], errors='coerce')
        kcat_values = df[column_name].dropna()

        total = len(df)
        matched = len(kcat_values)
        match_percent = matched / total * 100 if total else 0

        if not kcat_values.empty:
            min_exp = int(np.floor(np.log10(max(1e-6, kcat_values.min()))))
            max_exp = int(np.ceil(np.log10(kcat_values.max())))
            bins = np.logspace(min_exp, max_exp, num=40)

            fig, ax = plt.subplots(figsize=(8, 5))
            ax.hist(kcat_values, bins=bins, color="steelblue",
                    edgecolor="white", linewidth=0.7)

            ax.set_xscale("log")
            ax.set_xlim([10**min_exp / 1.5, 10**max_exp * 1.5])
            ax.xaxis.set_major_formatter(LogFormatter(10))

            ax.set_xlabel("kcat (s⁻¹)", fontsize=12)
            ax.set_ylabel("Count", fontsize=12)
            ax.set_title(f"{title} (n={matched}, {match_percent:.1f}%)", fontsize=13)

            return fig_to_base64(fig), matched, match_percent
        return "<p>No valid values available for plotting.</p>", 0, 0

    img_source, n_source, pct_source = plot_kcat_distribution(
        'kcat_source', "Experimental kcat distribution"
    )
    img_pred, n_pred, pct_pred = plot_kcat_distribution(
        'catapro_predicted_kcat_s', "Predicted kcat distribution"
    )

    # === 2. Difference boxplot ===
    if "kcat_source" in df.columns and "catapro_predicted_kcat_s" in df.columns:
        df["kcat_diff"] = df["catapro_predicted_kcat_s"] - df["kcat_source"]
        fig, ax = plt.subplots(figsize=(8, 5))
        df.boxplot(column="kcat_diff", by="matching_score", grid=False,
                   boxprops=dict(color="skyblue"), medianprops=dict(color="black"), ax=ax)
        ax.set_yscale("symlog")
        ax.set_ylabel("Difference (Predicted - Source kcat, s⁻¹)")
        ax.set_xlabel("Matching Score")
        ax.set_title("Difference of retrieved vs predicted kcat")
        plt.suptitle("")
        ax.axhline(0, color="red", linestyle="--", linewidth=1)
        img_diff = fig_to_base64(fig)
    else:
        img_diff = "<p>Required columns missing for difference plot.</p>"

    # === 3. Single segmented DB coverage bar ===
    db_counts = df["kcat_db"].fillna("Unknown").value_counts()
    total_db = db_counts.sum()

    # Improved color palette for better distinction and accessibility
    # Only 4 cases: brenda, sabio_rk, catapro, or None
    colors = {
    "brenda": "#3498db",      # Blue
    "sabio_rk": "#e67e22",    # Orange
    "catapro": "#2ecc71",     # Green
    "Unknown": "#7f8c8d"      # Gray (for None or unknown)
    }

    db_colors = {db: colors.get(db, "#7f8c8d") for db in db_counts.index}

    progress_segments = ""
    legend_items = ""
    for db, count in db_counts.items():
        percent = count / total_db * 100
        progress_segments += f"""
            <div class="progress-segment" style="width:{percent:.1f}%; background-color:{db_colors[db]};"
                 title="{db}: {percent:.1f}%"></div>
        """
        legend_items += f"""
            <span style="display:inline-block; margin-right:10px;">
                <span style="display:inline-block; width:15px; height:15px; background:{db_colors[db]}; margin-right:5px;"></span>
                {db} ({percent:.1f}%)
            </span>
        """

    progress_bar = f"""
        <div class="progress-multi">
            {progress_segments}
        </div>
        <div style="margin-top:10px;">{legend_items}</div>
    """

    # === 4. HTML layout with explanations ===
    html = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <title>WILDkCAT Report</title>
        <style>
            body {{
                font-family: 'Segoe UI', Tahoma, sans-serif;
                background-color: #f4f6f9;
                margin: 0;
                padding: 0;
                color: #333;
                line-height: 1.6;
            }}
            header {{
                background: linear-gradient(90deg, #2c3e50, #2980b9);
                color: #fff;
                padding: 20px;
                text-align: center;
            }}
            .container {{
                max-width: 1000px;
                margin: 20px auto;
                padding: 20px;
            }}
            .card {{
                background: #fff;
                border-radius: 10px;
                padding: 20px;
                margin-bottom: 20px;
                box-shadow: 0 2px 6px rgba(0,0,0,0.1);
            }}
            .card h2 {{
                color: #2980b9;
                margin-top: 0;
            }}
            .plot-container {{
                text-align: center;
                margin: 20px 0;
            }}
            .plot-container img {{
                max-width: 90%;
                border-radius: 6px;
                box-shadow: 0 2px 6px rgba(0,0,0,0.1);
            }}
            .progress-multi {{
                display: flex;
                width: 100%;
                height: 25px;
                border-radius: 12px;
                overflow: hidden;
                border: 1px solid #ccc;
            }}
            .progress-segment {{
                height: 100%;
            }}
        </style>
    </head>
    <body>
        <header>
            <h1>WILDkCAT Report</h1>
            <p>Generated on {generated_time}</p>
        </header>
        <div class="container">
            <div class="card">
                <h2>Introduction</h2>
                <p>
                    Lorem Ipsum
                </p>
            </div>

            <div class="card">
                <h2>Experimental k<sub>cat</sub> Distribution</h2>
                {img_source}
                <p>
                    Lorem Ipsum
                </p>
            </div>

            <div class="card">
                <h2>Predicted k<sub>cat</sub> Distribution</h2>
                {img_pred}
                <p>
                    Lorem Ipsum
                </p>
            </div>

            <div class="card">
                <h2>Difference Analysis</h2>
                {img_diff}
                <p>
                    Lorem Ipsum
                </p>
            </div>

            <div class="card">
                <h2>Coverage</h2>
                <p>
                    Lorem Ipsum
                </p>
                {progress_bar}
            </div>
        </div>

    <footer>WILDkCAT</footer>
    </body>
    </html>
    """

    os.makedirs("reports", exist_ok=True)
    report_path = "reports/general_report.html"
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(html)

    logging.info(f"HTML report saved to '{report_path}'")
    return report_path