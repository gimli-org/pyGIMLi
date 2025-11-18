import pygimli as pg

# log = logging.getLogger('pyGIMLi')

# logging.basicConfig(level=logging.DEBUG,
#                     format='%(asctime)s - %(name)s - ' + \
#                            '%(levelname)s - %(message)s',
#                     datefmt='%m/%d/%Y %H:%M:%S',
#                     #filename='example.log'
#                     )
pg.version()

# test pygimli log
pg.info("Start numeric log test." + str(pg.log(pg.Vector(1, 1.))))
pg.setVerbose(True)
pg.verbose("some verbose notes")
pg.warn("Start warning test.")


def testTraceback1():
    def testTraceback2():
        pg.error("Start error test.: int", 1, " vec", pg.Vector(2))
    testTraceback2()


testTraceback1()


@pg.v
def testVerboseDecorator1():
    pg.verbose('testVerboseDecorator1 should be seen even if verbose is false')


def testVerboseDecorator2():
    pg.verbose('testVerboseDecorator2 should be seen even if verbose is true')


pg.setVerbose(False)
testVerboseDecorator1()

pg.setVerbose(False)
testVerboseDecorator2()
pg.setVerbose(True)
testVerboseDecorator2()


@pg.d
def testDebugDecorator():
    pg.debug('testDebugDecorator should always be seen even if debug is false')


testDebugDecorator()

pg.d(testDebugDecorator())

pg.debug("debug 0")

pg.setDebug(1)
pg.debug("debug ON")
pg.setThreadCount(2)

pg.debug("debug with a trace", withTrace=True)


# should not printed out
pg.setDebug(0)
pg.debug("debug OFF")
pg.setThreadCount(2)

# test core log (should not be used outside the core)
pg.log(pg.core.Info, "core log ")
pg.log(pg.core.Warning, "core log ")
pg.log(pg.core.Error, "core log ")
try:
    pg.log(pg.core.Critical, "core log ")
except BaseException as e:
    print("cought core exception:", e)
# pg.exception("Exception")


def testMethod(**kwargs):
    pg.core.logger.warnNonEmptyArgs(kwargs)


testMethod(a=1, b='foo')
#

# teste colored output
print(pg._('Green', c='g'), pg._('Red', c='r'), pg._('Yellow', c='y'))
# print(pg._g('Green'), pg._g('Red'), pg._y('Yellow'))
print(pg._('more', 'then', 'one', c='6;30;42'))
