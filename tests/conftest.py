"""How Hypothesis draws its examples."""

from __future__ import annotations

import os

from hypothesis import settings

# On a developer's machine the property tests draw new examples on every run, which is how they find new faults.
# In CI the draw is fixed, so a red run means the code changed and not the dice.
settings.register_profile("ci", derandomize=True)
if os.environ.get("CI"):
    settings.load_profile("ci")
