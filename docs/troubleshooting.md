# Troubleshooting 

This section summarizes common issues and how to resolve them. 

## API Connection Issues

**Problem:** The SABIO-RK API is unavailable or unresponsive.

<details>
<summary>Example error</summary>

```bash
❯ python -m scripts.run_wildkcat_perso
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

