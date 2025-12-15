import requests
import logging 
from time import sleep
from functools import wraps


def retry_api(max_retries=4, backoff_factor=2):
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            delay = 2
            for attempt in range(max_retries):
                try:
                    return func(*args, **kwargs)

                except requests.exceptions.ConnectionError as e:
                    logging.warning(
                        f"Connection error ({e}), retry {attempt+1}/{max_retries}"
                    )

                except requests.exceptions.HTTPError as e:
                    status = e.response.status_code
                    if status in {502, 503, 504}:
                        logging.warning(
                            f"HTTP {status} error, retry {attempt+1}/{max_retries}"
                        )
                    else:
                        logging.error(f"HTTP error (no retry): {e}")
                        return None

                sleep(delay)
                delay *= backoff_factor

            logging.error(f"Request failed after {max_retries} retries.")
            return None

        return wrapper
    return decorator

def safe_requests_get(url, timeout=10):
    response = requests.get(url, timeout=timeout)
    response.raise_for_status()
    return response