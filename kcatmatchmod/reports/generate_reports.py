import logging
import pandas as pd
import datetime

# TODO: Generation of the reports of the different API are currently the same, it can be only one function that takes the DataFrame and the API name as parameters.


def report_extraction(model, df, transferred):
    """Generate a HTML report summarizing the kcat extraction results."""
    # Model statistics
    num_model_reactions = len(model.reactions)
    num_model_metabolites = len(model.metabolites)
    num_model_genes = len(model.genes)

    # kcat output statistics
    num_reactions = df['rxn'].nunique()
    num_ec_codes = df['ec_code'].nunique()
    num_ec_codes_transferred = transferred
    num_ec_codes_in_model = num_ec_codes + num_ec_codes_transferred
    num_kegg_rxn_ids = df['KEGG_rxn_id'].nunique()

    rxn_coverage = 100.0 * num_reactions / num_model_reactions if num_model_reactions else 0
    percent_ec_retrieved = 100.0 * num_ec_codes / (num_ec_codes + num_ec_codes_transferred) if (num_ec_codes + num_ec_codes_transferred) else 0

    # Timestamp
    generated_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    html = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Kcat Extraction Report</title>
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
            <h1>Kcat Extraction Report</h1>
            <p>Generated on {generated_time}</p>
        </header>

        <div class="container">
            <div class="card">
                <h2>Model Overview</h2>
                <div class="stats-grid">
                    <div class="stat-box">
                        <h3>{model.id}</h3>
                        <p>Model ID</p>
                    </div>
                    <div class="stat-box">
                        <h3>{num_model_reactions}</h3>
                        <p>Reactions</p>
                    </div>
                    <div class="stat-box">
                        <h3>{num_model_metabolites}</h3>
                        <p>Metabolites</p>
                    </div>
                    <div class="stat-box">
                        <h3>{num_model_genes}</h3>
                        <p>Genes</p>
                    </div>
                </div>
            </div>

            <div class="card">
                <h2>Kcat Extraction Statistics</h2>
                <table>
                    <tr>
                        <th>Metric</th>
                        <th>Value</th>
                        <th>Visualization</th>
                    </tr>
                    <tr>
                        <td>Reactions with EC info</td>
                        <td>{num_reactions} ({rxn_coverage:.1f}%)</td>
                        <td>
                            <div class="progress">
                                <div class="progress-bar" style="width:{rxn_coverage}%;"></div>
                            </div>
                        </td>
                    </tr>
                    <tr>
                        <td>EC codes found in KEGG</td>
                        <td>{num_ec_codes} (100%)</td>
                        <td>
                            <div class="progress">
                                <div class="progress-bar" style="width:100%;"></div>
                            </div>
                        </td>
                    </tr>
                    <tr>
                        <td>EC codes not found (transferred)</td>
                        <td>{num_ec_codes_transferred}</td>
                        <td>-</td>
                    </tr>
                    <tr>
                        <td>Unique KEGG reaction IDs</td>
                        <td>{num_kegg_rxn_ids}</td>
                        <td>-</td>
                    </tr>
                    <tr>
                        <td>Total rows in output</td>
                        <td>{len(df)}</td>
                        <td>-</td>
                    </tr>
                </table>
            </div>
        </div>

        <footer>
            Report generated automatically. &copy; {datetime.datetime.now().year}
        </footer>
    </body>
    </html>
    """

    # Save the report
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


def report_brenda(df): 
    """
    Generate an HTML report summarizing the kcat matching results from BRENDA.

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
    <head><title>BRENDA kcat Matching Report</title></head>
    <body>
        <h1>BRENDA kcat Matching Report</h1>
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
    report_path = "reports/brenda_report.html"
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(html)
    logging.info(f"HTML report saved to '{report_path}'")
