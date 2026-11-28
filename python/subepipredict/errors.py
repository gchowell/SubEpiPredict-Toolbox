"""Errors and warnings that carry the MATLAB message identifiers.

The MATLAB toolbox raises errors such as
``error('SubEpiPredict:InvalidModelFlag', ...)`` and warnings through
``subepidemicWarning.m``, which also records them so that a driver can list
them again at its end (``subepidemicWarningLog.m``). The Python port keeps the
identifier on the exception and on the warning (``exc.identifier``), and the
same log.
"""

from __future__ import annotations

import warnings


class SubEpiPredictError(Exception):
    """Error with a MATLAB-style identifier, e.g. ``SubEpiPredict:InvalidModelFlag``."""

    def __init__(self, identifier, message=""):
        self.identifier = identifier
        self.message = message
        super().__init__("%s: %s" % (identifier, message) if message else identifier)


class SubEpiPredictWarning(UserWarning):
    """Warning with a MATLAB-style identifier (``w.identifier``)."""

    def __init__(self, identifier, message=""):
        self.identifier = identifier
        self.message = message
        super().__init__("%s: %s" % (identifier, message) if message else identifier)


# The warnings of one run: list of dicts with the keys id, message and count.
_ENTRIES = []


def _switched_off(warning):
    """True when the warnings filters would ignore this warning."""
    text = str(warning)
    for action, message, category, _module, _lineno in warnings.filters:
        if not issubclass(type(warning), category):
            continue
        if message is not None and not message.match(text):
            continue
        return action == "ignore"
    return False


def subepidemicWarning(identifier, messageFormat, *args, stacklevel=2):
    """Raise a warning of the toolbox and remember it for the summary.

    Port of ``subepidemicWarning.m``. The message is ``messageFormat % args``.
    A warning that the ``warnings`` filters ignore is neither shown nor
    recorded.
    """
    message = messageFormat % args if args else messageFormat
    w = SubEpiPredictWarning(identifier, message)
    if not _switched_off(w):
        subepidemicWarningLog("add", identifier, message)
    warnings.warn(w, stacklevel=stacklevel + 1)


def subepidemicWarningLog(action, *args, file=None):
    """Keep the toolbox's warnings of one run and print them at its end.

    Port of ``subepidemicWarningLog.m``.

    ``subepidemicWarningLog('start')``           empties the list
    ``subepidemicWarningLog('add', id, message)`` called by :func:`subepidemicWarning`
    ``subepidemicWarningLog('summary', label)``  prints the warnings raised since
                                                 'start', one line per kind
    ``subepidemicWarningLog('get')``             returns the list (dicts with the
                                                 keys ``id``, ``message``, ``count``)
    """
    global _ENTRIES
    if action == "start":
        _ENTRIES = []
    elif action == "add":
        identifier, message = args[0], args[1]
        for entry in _ENTRIES:
            if entry["id"] == identifier:
                entry["count"] += 1
                break
        else:
            _ENTRIES.append({"id": identifier, "message": message, "count": 1})
    elif action == "summary":
        label = args[0] if args else "SubEpiPredict"
        if _ENTRIES:
            print("\n%s finished with %d kind(s) of warning:" % (label, len(_ENTRIES)), file=file)
            for entry in _ENTRIES:
                if entry["count"] == 1:
                    print("  %s: %s" % (entry["id"], entry["message"]), file=file)
                else:
                    print("  %s (%d times; first:) %s" % (entry["id"], entry["count"], entry["message"]), file=file)
    elif action == "get":
        return [dict(entry) for entry in _ENTRIES]
    else:
        raise SubEpiPredictError("SubEpiPredict:WarningLog", 'Unknown action "%s".' % action)
    return None
