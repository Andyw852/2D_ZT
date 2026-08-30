"""Phase A - Step 0: 记录 Python 环境信息。"""
import sys
import platform
from pathlib import Path
from datetime import datetime, timezone

import jarvis
import numpy
import pandas
import requests

root = Path(__file__).resolve().parents[1]
reports = root / "reports"
reports.mkdir(exist_ok=True)

now = datetime.now(timezone.utc).isoformat()
lines = [
    f"execution_date={now}",
    f"python_version={platform.python_version()}",
    f"python_executable={sys.executable}",
    f"platform={platform.platform()}",
    f"jarvis_tools_version={getattr(jarvis, '__version__', 'unknown')}",
    f"numpy_version={numpy.__version__}",
    f"pandas_version={pandas.__version__}",
    f"requests_version={requests.__version__}",
]
(reports / "python_version.txt").write_text(
    "\n".join(lines) + "\n", encoding="utf-8"
)
print("\n".join(lines))
