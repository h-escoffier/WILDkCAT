import logging
import os
import datetime
import logging
import base64
import io
import numpy as np
import matplotlib.pyplot as plt


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
            <h1>Kcat Extraction Report</h1>
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

            <!-- Kcat Extraction Table -->
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
                        <td>Number of kcat values dropped due to inconsistent EC codes</td>
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

        <footer>
            KcatMetaMod
        </footer>
    </body>
    </html>
    """

    # Save report
    os.makedirs("reports", exist_ok=True)
    report_path = "reports/kcat_summary.html"
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(html)
    logging.info(f"HTML report saved to '{report_path}'")



def report_api(df, api_name):
    """
    Generate a styled HTML report summarizing the kcat matching results from SABIO-RK/BRENDA,
    including kcat value distribution and matching score repartition.

    Parameters:
        df (pd.DataFrame): DataFrame containing columns ['kcat', 'matching_score'].
        api_name (str): Name of the API (e.g., 'SABIO-RK', 'BRENDA').
    """

    # Gather statistics
    total = len(df)
    matched = df['kcat'].notna().sum()
    match_percent = matched / total * 100 if total > 0 else 0

    # Ensure all scores 1-5 and 10 are present in the index (even if 0)
    all_scores = [1, 2, 3, 4, 5, 10]
    score_counts = df['matching_score'].value_counts().reindex(all_scores, fill_value=0)
    score_percent = (score_counts / total * 100).round(2)

    if api_name.lower() == 'sabio_rk':
        api = 'SABIO-RK'
    elif api_name.lower() == 'brenda':
        api = 'BRENDA'
    else:
        api = api_name
        logging.error(f"Unknown API: {api_name}")

    # kcat value stats
    kcat_values = df['kcat'].dropna()

    kcat_hist_base64 = ""
    if not kcat_values.empty:
        # Use log10 bins from 10^-1 to 10^2 (or adapt to data range)
        min_exp = max(-1, int(np.floor(np.log10(kcat_values[kcat_values > 0].min())))) if (kcat_values > 0).any() else -1
        max_exp = min(2, int(np.ceil(np.log10(kcat_values.max())))) if (kcat_values > 0).any() else 1
        bins = np.logspace(min_exp, max_exp, num=40)
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.hist(kcat_values, bins=bins, color="#2980b9", alpha=0.85, edgecolor='white')
        ax.set_xscale('log')
        ax.set_xlabel("kcat", fontsize=14)
        ax.set_ylabel("Count", fontsize=14)
        ax.tick_params(axis='both', labelsize=12)
        ax.set_xlim([10**min_exp / 1.5, 10**max_exp * 1.5])
        ax.set_xticks([10**i for i in range(min_exp, max_exp+1)])
        ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
        plt.tight_layout()
        buf = io.BytesIO()
        plt.savefig(buf, format='png', bbox_inches='tight')
        plt.close(fig)
        kcat_hist_base64 = base64.b64encode(buf.getvalue()).decode('utf-8')

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
            .bar-1 {{ background-color: #2980b9; }}  /* Perfect match */
            .bar-2 {{ background-color: #16a085; }}  /* Relax enzyme */
            .bar-3 {{ background-color: #f39c12; }}  /* Relax pH & temperature */
            .bar-4 {{ background-color: #8e44ad; }}  /* Relax organism */
            .bar-5 {{ background-color: #c0392b; }}  /* Relax substrate */
            .bar-10 {{ background-color: #ecf0f1; color: #333; }}  /* Not found (light grey) */
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
            .img-section {{
                display: flex;
                flex-wrap: wrap;
                gap: 30px;
                justify-content: center;
                align-items: flex-start;
                margin-top: 20px;
            }}
            .img-card {{
                background: #f9fafc;
                border-radius: 8px;
                padding: 15px;
                border: 1px solid #e2e2e2;
                text-align: center;
                margin-bottom: 10px;
            }}
            .img-card img {{
                max-width: 350px;
                display: block;
                margin: 0 auto 10px auto;
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

    # Color and label maps for scores 1-5 and 10
    color_map = {
        1: "bar-1",  # Perfect match
        2: "bar-2",  # Relax enzyme
        3: "bar-3",  # Relax pH & temperature
        4: "bar-4",  # Relax organism
        5: "bar-5",  # Relax substrate
        10: "bar-10" # Not found (light grey)
    }
    label_map = {
        1: "Perfect match",
        2: "Relax enzyme",
        3: "Relax pH & temperature",
        4: "Relax organism",
        5: "Relax substrate",
        10: "Not found"
    }

    # Add each segment for scores 1-5 and 10
    for score in all_scores:
        percent = score_percent.get(score, 0)
        bar_class = color_map.get(score, "bar-10")
        label = label_map.get(score, f"Score {score}")
        if percent > 0:
            html += f"""
                    <div class="progress-bar {bar_class}" style="width:{percent}%" title="{label}: {percent:.2f}%"></div>
            """

    html += """
                </div>
                <div class="legend">
    """

    for score in all_scores:
        bar_class = color_map.get(score, "bar-10")
        label = label_map.get(score, f"Score {score}")
        html += f"""
                    <div class="legend-item">
                        <div class="legend-color {bar_class}"></div> {label}
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
    for score in all_scores:
        count = score_counts[score]
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
    """

    html += """
            <div class="card">
                <h2>Distribution of kcat values</h2>
                <div class="img-section">
    """
    if kcat_hist_base64:
        html += f"""
        <img src="data:image/png;base64,{kcat_hist_base64}" alt="kcat Distribution">
    """
    
    html += """
                </div>
            </div>
    """

    html += """
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
            KcatMetaMod
        </footer>
    </body>
    </html>
    """

    # Save to file
    os.makedirs("reports", exist_ok=True)
    report_path = f"reports/{api_name}_report.html"
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(html)

    logging.info(f"HTML report saved to '{report_path}'")
