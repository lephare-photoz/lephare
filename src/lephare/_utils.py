import sys

__all__ = [
    "continueClass",
]

# code copied from https://github.com/lsst/utils/blob/main/python/lsst/utils/wrappers.py
INTRINSIC_SPECIAL_ATTRIBUTES = frozenset(
    (
        "__qualname__",
        "__module__",
        "__metaclass__",
        "__dict__",
        "__weakref__",
        "__class__",
        "__subclasshook__",
        "__name__",
        "__doc__",
    )
)


# code copied from https://github.com/lsst/utils/blob/main/python/lsst/utils/wrappers.py
def isAttributeSafeToTransfer(name, value):  # noqa: N802
    """Return True if an attribute is safe to monkeypatch-transfer to another
    class.
    This rejects special methods that are defined automatically for all
    classes, leaving only those explicitly defined in a class decorated by
    `continueClass` or registered with an instance of `TemplateMeta`.
    """
    if name.startswith("__") and (
        value is getattr(object, name, None) or name in INTRINSIC_SPECIAL_ATTRIBUTES
    ):
        return False
    return True


# code copied from https://github.com/lsst/utils/blob/main/python/lsst/utils/wrappers.py
def continueClass(cls):  # noqa: N802
    orig = getattr(sys.modules[cls.__module__], cls.__name__)
    for name in dir(cls):
        # Common descriptors like classmethod and staticmethod can only be
        # accessed without invoking their magic if we use __dict__; if we use
        # getattr on those we'll get e.g. a bound method instance on the dummy
        # class rather than a classmethod instance we can put on the target
        # class.
        attr = cls.__dict__.get(name, None) or getattr(cls, name)
        if isAttributeSafeToTransfer(name, attr):
            setattr(orig, name, attr)
    return orig
