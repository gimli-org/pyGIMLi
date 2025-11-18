# -*- coding: utf-8 -*-
"""pyGIMLi - Stuff for logging

    use with
    pg.info()
    pg.warn()
    pg.debug()
    pg.error() # non critical
    pg.critical() # raises exception

"""
import logging
import inspect
import traceback

from . core import pgcore

__ANSICOLORS__ = {
    'r': '\033[0;31;49m',  # normal, #FG red; #BG black
    'g': '\033[0;32;49m',
    'y': '\033[0;33;49m',
    'b': '\033[0;34;49m',
    'm': '\033[0;35;49m',
    'DARK_GREY': '\033[1;30m',
    'NC': '\033[0m',  # No Color
}


def _msg(*args):
    msg = ''
    for i, arg in enumerate(args):
        msg += str(arg)
        if i < len(args)-1:
            msg += ' '
    return msg


def _(*args, c=None):
    # will probably only for linux or any msys like shell
    if c is None:
        return _msg(*args)
    elif '\033' in c:
        return c + _msg(*args) + __ANSICOLORS__['NC']
    try:
        return __ANSICOLORS__[c] + _msg(*args) + __ANSICOLORS__['NC']
    except Exception:
        return '\033[' + c + 'm' + _msg(*args) + __ANSICOLORS__['NC']


def _get_class_from_frame(fr):
    args, _, _, value_dict = inspect.getargvalues(fr)
    if len(args) and args[0] == 'self':
        instance = value_dict.get('self', None)
        if instance is not None:
            return getattr(instance, '__class__', None)
    return None


def whereAmI(nr=3):
    nr = min(len(inspect.stack())-1, nr)
    clsName = _get_class_from_frame(inspect.stack()[nr][0])
    method = inspect.stack()[nr].function
    fileN = inspect.stack()[nr].filename.split('/')[-1]
    line = inspect.stack()[nr].lineno
    # print(inspect.stack()[nr])
    # print('{0}:{1}'.format(fileN, line))
    return str(clsName) + '.' + method + '({0}:{1})'.format(fileN, line)


def p(*args, c='y'):
    deprecated(hint='use _d instead')
    print(_(*args, c=c))


def _g(*args):
    _d(*args, c='g')


def _y(*args):
    _d(*args, c='y')


def _r(*args):
    _d(*args, c='r')


def _b(*args):
    _d(*args, c='b')


def _d(*args, c='y'):
    """Simplistic colored debug msg"""
    print(_(whereAmI(), ':', *args, c=c))


class ColorFormatter(logging.Formatter):
    def __init__(self,
                 fmt='%(asctime)s - %(name)s - %(levelname)s - %(message)s'):
        logging.Formatter.__init__(self, fmt, "%d/%m/%y - %H:%M:%S")
        self._formatter = {
            logging.INFO: logging.PercentStyle(
                '%(asctime)s - %(name)s - ' + _('%(levelname)s', c='g') +
                ' - %(message)s'),
            logging.ERROR: logging.PercentStyle(
                '%(asctime)s - %(name)s - ' + _('%(levelname)s', c='r') +
                ' - %(message)s'),
            logging.WARNING: logging.PercentStyle(
                '%(asctime)s - %(name)s - ' + _('%(levelname)s', c='y') +
                ' - %(message)s'),
            logging.CRITICAL: logging.PercentStyle(
                '%(asctime)s - %(name)s - ' +
                _('%(levelname)s', c='5;31;1;49') + ' - %(message)s'),
            logging.DEBUG: logging.PercentStyle(
                '%(asctime)s - %(name)s - ' + _('%(levelname)s', c='m') +
                ' - %(message)s'),
            'DEFAULT': logging.PercentStyle(
                '%(asctime)s - %(name)s - %(levelname)s - %(message)s'),
        }

    def format(self, record):
        self._style = self._formatter.get(record.levelno,
                                          self._formatter['DEFAULT'])
        return logging.Formatter.format(self, record)


logger = logging.getLogger('pyGIMLi')
streamHandler = logging.StreamHandler()
streamHandler.setFormatter(ColorFormatter())
logger.root.addHandler(streamHandler)


def setLogLevel(level):
    """Shortcut to change the current log level"""
    logger.setLevel(level)


def addLogLevel(value, name):
    """Add a new log level to the :mod:`logging` module.

    Parameters
    ----------
    value: int
        log level number.
    name: str
        Name for the log level.
    """
    logging.addLevelName(value, name)
    setattr(logging, name, value)


# CRITICAL = 50
# ERROR = 40
# WARNING = 30
# INFO = 20
VERBOSE = 15
# DEBUG = 10
# NOTSET = 0
addLogLevel(VERBOSE, 'VERBOSE')


def __logVerbose(msg, *args, **kwargs):
    if logger.isEnabledFor(VERBOSE):
        logger._log(VERBOSE, msg, args, **kwargs)


logger.verbose = __logVerbose

__verbose_level__ = 0


class VerboseScope(object):
    def __init__(self, verb):
        self._oldState = __verbose_level__
        self._state = verb
        setVerbose(self._state)

    def __del__(self):
        setVerbose(self._oldState)

    def __bool__(self):
        if self._state == 1:
            return True
        else:
            return False


def setVerbose(v):
    level = logging.INFO
    if v:
        __verbose_level__ = 1
        level = logging.VERBOSE
    else:
        __verbose_level__ = 0
        level = logging.INFO

    logger.setLevel(level)


def v(funct):
    """Decorator to enable verbose messages for the scope of a function.

    Examples
    --------
    >>> import pygimli as pg
    >>> @pg.v                                                  # doctest: +SKIP
    >>> def foo():                                             # doctest: +SKIP
    ...     pg.verbose('foo')                                  # doctest: +SKIP
    >>> foo()                                                  # doctest: +SKIP
    >>> def bar(d):
    ...     pg.verbose('bar', d)
    >>> bar('verbose should be off')
    >>> pg.setVerbose(1)
    >>> bar('verbose should be on (1)')
    >>> pg.setVerbose(0)
    >>> pg.v(bar)('verbose should be on (2)')
    """
    def wrapper(*args, **kwargs):
        o = logger.level
        logger.setLevel(logging.VERBOSE)
        rv = funct(*args, **kwargs)
        logger.setLevel(o)
        return rv
    return wrapper


def d(funct):
    """Decorator to enable debug messages for the scope of a function.

    Examples
    --------
    >>> import pygimli as pg
    >>> @pg.d                                                  # doctest: +SKIP
    >>> def foo():                                             # doctest: +SKIP
    ...     pg.debug('foo')                                    # doctest: +SKIP
    >>> foo()                                                  # doctest: +SKIP
    >>> def bar(d):
    ...     pg.debug('bar', d)
    >>> bar('debug should be off')
    >>> pg.setDebug(1)
    >>> bar('debug should be on (1)')
    >>> pg.setDebug(0)
    >>> pg.d(bar)('debug should be on (2)')
    """
    def wrapper(*args, **kwargs):
        o = logger.level
        logger.setLevel(logging.DEBUG)
        rv = funct(*args, **kwargs)
        logger.setLevel(o)
        return rv
    return wrapper


def verbose():
    return __verbose_level__


def setDebug(d):
    level = logging.INFO
    if d:
        pgcore.setDebug(True)
        level = logging.DEBUG
    else:
        pgcore.setDebug(False)
        level = logging.INFO

    logger.setLevel(level)
    logging.getLogger('Core').setLevel(level)
    fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=level, format=fmt, datefmt='%m/%d/%Y %H:%M:%S')
    # filename='pygimli.log'


def info(*args):
    logger.info(_msg(*args))


def warn(*args):
    logger.warning(_msg(*args))


def error(*args):
    logger.error(whereAmI(nr=2) + "\n" + _msg(*args))


def debug(*args, withTrace=False):
    """
    Parameters
    ----------
    """
    if withTrace:
        traceback.print_exc()
    logger.debug(_msg(*args))


def verbose(*args):   # isn't this a refinition of line 253?
    logger.verbose(_msg(*args))


def critical(*args):
    logger.critical(whereAmI(nr=2) + "\n" + _msg(*args))
    raise Exception(_msg(*args))


def deprecated(msg='', hint=''):
    logger.warning("Deprecated code usage at:")
    logger.warning(whereAmI() + "\n" + msg + " " + hint)


# def renamed(newFunc, removed=''):
#     """Rename the current function into newFunc.
#     Give remove date or version nummer.
#     """
#     logger.warning(whereAmI() + ' is renamed into ' + newFunc.__name__ +
#                    ' and will be removed in: ' + removed)
#     return newFunc(**inspect.stack()[1].frame.f_locals)


def renameKwarg(old, new, kwargs, ver=''):
    if old in kwargs:
        logger.warning("Keyword argument changed from '" + old +
                       "' to '" + new + "' and will be removed in v " + ver)
        kwargs[new] = kwargs.pop(old)


def warnNonEmptyArgs(kwargs):
    if len(kwargs) > 0:
        logger.warning(whereAmI() +
                       "Unrecognized keyword arguments for method:" +
                       _msg(kwargs))
