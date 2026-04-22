# Troubleshooting 

This section summarizes common issues and how to resolve them. 

## API Connection Issues

**Problem:** The SABIO-RK API is unavailable or unresponsive.

<details>
<summary>Example error</summary>

```bash
Retrieving kcat values:   0%|                                                                                                                                                                                                                          | 0/6072 [00:01<?, ?it/s]
Traceback (most recent call last):
  File "<frozen runpy>", line 198, in _run_module_as_main
  File "<frozen runpy>", line 88, in _run_code
  File "/Users/hugues_esc/VSCodeProjects/KcatMatchMod/scripts/run_wildkcat_perso.py", line 11, in <module>
    run_retrieval(
  File "/Users/hugues_esc/VSCodeProjects/KcatMatchMod/wildkcat/processing/retrieve_kcat.py", line 260, in run_retrieval
    best_match, matching_score = extract_kcat(kcat_dict, general_criteria, database=database)
                                 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/Users/hugues_esc/VSCodeProjects/KcatMatchMod/wildkcat/processing/retrieve_kcat.py", line 118, in extract_kcat
    api_output = get_turnover_number(kcat_dict['ec_code'], database)
                 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/Users/hugues_esc/VSCodeProjects/KcatMatchMod/wildkcat/processing/retrieve_kcat.py", line 42, in get_turnover_number
    df_sabio = get_turnover_number_sabio(ec_code)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/Users/hugues_esc/VSCodeProjects/KcatMatchMod/wildkcat/api/sabio_rk_api.py", line 30, in get_turnover_number_sabio
    request.raise_for_status()
  File "/opt/anaconda3/envs/wildkcat-env/lib/python3.11/site-packages/requests/models.py", line 1026, in raise_for_status
    raise HTTPError(http_error_msg, response=self)
requests.exceptions.HTTPError: 404 Client Error:  for url: https://sabiork.h-its.org/sabioRestWebServices/searchKineticLaws/entryIDs?format=txt&q=Parametertype%3A%22kcat%22+AND+ECNumber%3A%222.7.1.48%22
```

</details>

**Cause:** The Sabio-RK service may be temporarily down or unreachable.

**Solution:** 

1. Verify that the Sabio-RK website is accessible: [https://sabiork.h-its.org](https://sabiork.h-its.org)
2. If the site is online, try rerunning your script after a few minutes.
3. Otherwise, rerun the script with the `--database=brenda` option or wait until the Sabio-RK service is back online.

--- 

## `.env` Variables not loaded

**Problem:** Environment variables defined in the .env file are not being read by WILDkCAT, causing it to fail on startup or during execution.

<details>
<summary>Example error</summary>

```bash
Retrieving kcat values:   0%|                                                                                                                                                                                  | 0/226 [00:00<?, ?it/s]
Traceback (most recent call last):
  File "<frozen runpy>", line 198, in _run_module_as_main
  File "<frozen runpy>", line 88, in _run_code
  File "/Users/hugues_esc/VSCodeProjects/KcatMatchMod/scripts/run_wildkcat_perso.py", line 12, in <module>
    run_retrieval(
  File "/Users/hugues_esc/VSCodeProjects/KcatMatchMod/wildkcat/processing/retrieve_kcat.py", line 302, in run_retrieval
    best_match, penalty_score = extract_kcat(kcat_dict, general_criteria, database=database)
                                ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/Users/hugues_esc/VSCodeProjects/KcatMatchMod/wildkcat/processing/retrieve_kcat.py", line 119, in extract_kcat
    api_output = get_turnover_number(kcat_dict['ec_code'], database)
                 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/Users/hugues_esc/VSCodeProjects/KcatMatchMod/wildkcat/processing/retrieve_kcat.py", line 41, in get_turnover_number
    df_brenda = get_turnover_number_brenda(ec_code)
                ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/Users/hugues_esc/VSCodeProjects/KcatMatchMod/wildkcat/api/brenda_api.py", line 82, in get_turnover_number_brenda
    email, hashed_password = get_brenda_credentials()
                             ^^^^^^^^^^^^^^^^^^^^^^^^
  File "/Users/hugues_esc/VSCodeProjects/KcatMatchMod/wildkcat/api/brenda_api.py", line 62, in get_brenda_credentials
    raise ValueError("BRENDA_EMAIL and BRENDA_PASSWORD environment variables must be set.")
ValueError: BRENDA_EMAIL and BRENDA_PASSWORD environment variables must be set.
```

</details>

**Cause:** On Windows, Conda environments do not automatically load `.env` files. As a result, variables like `ENTREZ_EMAIL`, `BRENDA_EMAIL`, and `BRENDA_PASSWORD` are never injected into the environment, even if the `.env` file exists and is correctly formatted.

**Solution:** Set the variables directly in your Conda environment so they are loaded automatically on activation:

```bash
conda activate env-name
conda env config vars set ENTREZ_EMAIL=your_email@example.com
conda env config vars set BRENDA_EMAIL=your_email@example.com
conda env config vars set BRENDA_PASSWORD=your_password

# Re-activate the environment to apply the changes
conda activate env-name
```

--- 

