import requests
import logging 
from time import sleep
from functools import wraps


def retry_api(max_retries=3, backoff_factor=1.5):
    """Decorator to retry API calls with exponential backoff."""
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            delay = 1
            for attempt in range(max_retries):
                try:
                    return func(*args, **kwargs)
                except (requests.RequestException, ValueError) as e:
                    logging.warning(f"API call failed ({e}), retry {attempt+1}/{max_retries}...")
                    sleep(delay)
                    delay *= backoff_factor
            logging.error(f"API call failed after {max_retries} retries.")
            return None
        return wrapper
    return decorator


def safe_requests_get(url, timeout=10):
    response = requests.get(url, timeout=timeout)
    response.raise_for_status()
    return response