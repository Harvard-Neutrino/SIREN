"""Copy Python wrapper state without invoking archive guards."""

from copy import deepcopy
from types import MemberDescriptorType


def copy_state(source, memo=None):
    cls = type(source)
    result = object.__new__(cls)
    if memo is not None:
        memo[id(source)] = result
        result.__dict__ = deepcopy(source.__dict__, memo)
    else:
        result.__dict__ = source.__dict__.copy()
    for base in cls.__mro__:
        for member in vars(base).values():
            if isinstance(member, MemberDescriptorType):
                try:
                    value = member.__get__(source, cls)
                except AttributeError:
                    continue
                member.__set__(result, value if memo is None else deepcopy(value, memo))
    return result
