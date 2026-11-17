"""Minimal boost::python binding check used for debugging."""

print("#" * 60, "bp")

import _bpt_  # noqa: E402  deliberately after the banner

print(dir(_bpt_))  # shows the bound methods
print("1.0 =", _bpt_.t1())  # expecting 1 , works
print("2 =", _bpt_.t2())  # exception 2 , segfault
