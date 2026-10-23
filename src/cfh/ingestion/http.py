"""A ``requests`` session that retries transient network failures.

The cBioPortal fetchers already retry retryable HTTP status codes in their
own loops, but a timeout or a failed connection raises out of
``session.post`` on the first attempt. A cohort scan makes thousands of
calls to cBioPortal and Genome Nexus, so one such blip used to drop a gene
from the FDR family (#116). The adapter here retries connection errors and
read timeouts for every method; status codes are left to the callers so
retries do not multiply.
"""

from __future__ import annotations

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

DEFAULT_NETWORK_RETRIES = 6
DEFAULT_BACKOFF_FACTOR = 2.0


def retrying_session(
    *,
    retries: int = DEFAULT_NETWORK_RETRIES,
    backoff_factor: float = DEFAULT_BACKOFF_FACTOR,
) -> requests.Session:
    """Session retrying connection errors and read timeouts with backoff.

    Every endpoint the pipeline calls is a read, including the ``/fetch``
    POSTs, so retrying any method is safe.
    """
    retry = Retry(
        total=None,
        connect=retries,
        read=retries,
        other=retries,
        status=0,
        allowed_methods=None,
        backoff_factor=backoff_factor,
        backoff_max=60,
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry)
    session = requests.Session()
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session
