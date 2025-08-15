import os
import io
import base64
import logging
import datetime
import pandas as pd
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
    Generate a styled HTML report summarizing the kcat matching results,
    including kcat value distribution and matching score repartition.

    Parameters:
        df (pd.DataFrame): Must contain ['kcat', 'matching_score', ...].
        api_name (str): Name of the API (e.g., 'SABIO-RK', 'BRENDA').
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
    distinct_colors = [
        "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
        "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
        "#8CD17D", "#B6992D", "#499894", "#D37295", "#FABFD2",
        "#B2B2B2", "#5F9ED1", "#FFBE7D"
    ]
    def score_color(score):
        idx = present_scores.index(score)
        return distinct_colors[idx % len(distinct_colors)]

    generated_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Histogram with stacked bars for scores
    kcat_hist_base64 = ""
    if not kcat_values.empty:
        min_exp = max(-1, int(np.floor(np.log10(max(1e-6, kcat_values.min())))))
        max_exp = min(2, int(np.ceil(np.log10(kcat_values.max()))))
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
        ax.set_xticks([10**i for i in range(min_exp, max_exp+1)])
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
        <title>{api_name} kcat Matching Report</title>
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
            <h1>{api_name} kcat Matching Report</h1>
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
            <h2>Distribution of kcat values (Stacked by Matching Score)</h2>
            <div class="img-section">
    """
    if kcat_hist_base64:
        html += f'<img src="data:image/png;base64,{kcat_hist_base64}" alt="kcat Distribution">'
    html += """
            </div>
        </div>
    """

    # Metadata section
    html += """
            <div class="card">
                <h2>Matching Score Meaning</h2>
                <p>Matching score ranges from 0 (best match) to 15 (no match).</p>
            </div>
        </div>

        <footer>KcatMetaMod</footer>
    </body>
    </html>
    """

    # Save HTML
    os.makedirs("reports", exist_ok=True)
    report_path = f"reports/{api_name}_report.html"
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(html)

    logging.info(f"HTML report saved to '{report_path}'")


def report_catapro_part1(df, report_statistics):
    pass 


def report_catapro_part2(df, report_statistics):
    pass


def report_full(): 
    pass 
