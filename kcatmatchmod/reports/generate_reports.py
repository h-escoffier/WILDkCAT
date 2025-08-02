import logging
import os
import pandas as pd
import datetime


# TODO: The report_api function is not working well, problems with the precentage and the legend due to some missing values 


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
                        <td>{num_ec_codes} ({percent_ec_retrieved:.1f}%)</td>
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
        </div>

        <footer>
            KcatMetaMod - Report generated automatically. 
        </footer>
    </body>
    </html>
    """

    # Save the report
    report_path = "reports/kcat_summary.html"
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(html)
    logging.info(f"HTML report saved to '{report_path}'")


def report_api(df, api_name):
    """
    Generate a styled HTML report summarizing the kcat matching results from SABIO-RK.
    
    Parameters:
        df (pd.DataFrame): DataFrame containing columns ['kcat', 'matching_score'].
        api_name (str): Name of the API (e.g., 'SABIO-RK', 'BRENDA').
    """
    # Gather statistics
    total = len(df)
    matched = df['kcat'].notna().sum()
    match_percent = matched / total * 100 if total > 0 else 0
    score_counts = df['matching_score'].value_counts().sort_index()
    score_percent = (score_counts / total * 100).round(2)

    if api_name == 'sabio_rk':
        api = 'SABIO-RK'
    elif api_name == 'brenda':
        api = 'BRENDA'
    else:
        logging.error(f"Unknown API: {api_name}")

    # Current time
    generated_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Generate HTML content
    html = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>{api} kcat Matching Report</title>
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
            /* Stacked progress bar */
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
            /* Color palette */
            .bar-0 {{ background-color: #2980b9; }}  /* Perfect match */
            .bar-1 {{ background-color: #16a085; }}  /* Relax enzyme */
            .bar-2 {{ background-color: #f39c12; }}  /* Relax pH & temperature */
            .bar-3 {{ background-color: #8e44ad; }}  /* Relax organism */
            .bar-4 {{ background-color: #c0392b; }}  /* Relax substrate */
            .bar-5 {{ background-color: #ecf0f1; color: #333; }}  /* Not found (light grey) */
            /* Legend styles */
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
            }}
            ul {{
                padding-left: 20px;
                margin-top: 10px;
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
            <h1>{api} kcat Matching Report</h1>
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
                        <p>Matched kcat ({match_percent:.2f}%)</p>
                    </div>
                </div>
            </div>

            <div class="card">
                <h2>Matching Score Distribution</h2>
                <div class="progress-stacked">
    """

    # Add legend (same order as algorithm explanation)
    legend_labels = [
        "Perfect match",
        "Relax enzyme",
        "Relax pH & temperature",
        "Relax organism",
        "Relax substrate",
        "Not found"
    ]

    # Add each segment
    for i, label in enumerate(legend_labels):
        percent = score_percent.get(i, 0)  # use numeric index
        html += f"""
                    <div class="progress-bar bar-{i}" style="width:{percent}%" title="{label}: {percent:.2f}%"></div>
        """


    html += """
                </div>
                <div class="legend">
    """

    for i, label in enumerate(legend_labels):
        html += f"""
                    <div class="legend-item">
                        <div class="legend-color bar-{i}"></div> {label}
                    </div>
        """

    html += """
                </div>

                <table>
                    <tr>
                        <th>Score</th>
                        <th>Count</th>
                        <th>Percent</th>
                    </tr>
    """

    # Table rows
    for score, count in score_counts.items():
        percent = score_percent[score]
        html += f"""
                    <tr>
                        <td>{score}</td>
                        <td>{count}</td>
                        <td>{percent:.2f}%</td>
                    </tr>
        """

    html += """
                </table>
            </div>

            <div class="card">
                <h2>Hierarchical Matching Strategy</h2>
                <p>The algorithm attempts to find the best match in the following order:</p>
                <ul>
                    <li>1. Perfect match : EC - Enzyme - Temperature/pH - Organism - Substrate </li>
                    <li>2. Relax enzyme : EC - Temperature/pH - Organism - Substrate </li>
                    <li>3. Relax pH and temperature : EC - Organism - Substrate </li>
                    <li>4. Relax organism : EC - Substrate </li>
                    <li>5. Relax substrate : EC</li>
                    <li>10. Not found</li>
                </ul>
            </div>
        </div>

        <footer>
            KcatMetaMod - Report generated automatically.
        </footer>
    </body>
    </html>
    """


    # Save to file
    os.makedirs("reports", exist_ok=True)
    report_path = "reports/{api}_report.html".format(api=api_name)
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(html)

    logging.info(f"HTML report saved to '{report_path}'")
