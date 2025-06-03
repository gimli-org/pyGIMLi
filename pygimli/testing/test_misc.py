#!/usr/bin/env python
# -*- coding: utf-8 -*-

# write a correct test!
import unittest

import pygimli as pg
import numpy as np

class TestMisc(unittest.TestCase):

    def test_Trans(self):
        """
        """
        f = pg.trans.Trans()
        x = pg.Vector(3, 1.0)

        np.testing.assert_array_equal(f(x), x)
        np.testing.assert_array_equal(f.inv(x), x)
        np.testing.assert_array_equal(f.inv(f(x)), x)
        self.assertEqual(f.trans(1.0), 1.0)
        self.assertEqual(f(1.0), 1.0)
        self.assertEqual(f.inv(1.0), 1.0)

        f = pg.trans.TransLin(factor=2., offset=4.)
        np.testing.assert_array_equal(f(x), x*2. + 4.)
        np.testing.assert_array_equal(f.trans(x), x*2. + 4.)
        np.testing.assert_array_equal(f.inv(f(x)), x)
        self.assertEqual(f(1.0), 6.0)
        self.assertEqual(f.trans(1.0), 6.0)
        self.assertEqual(f.inv(6.0), 1.0)
        self.assertEqual(f.invTrans(6.0), 1.0)

        f = pg.trans.TransLogLU(lowerbound=0, upperbound=10)
        # print(f.update([1.], [100.]))
        np.testing.assert_array_equal(f.update([1.], [100.]), [10.0])
        # print(f.update([1.], [1000.]))
        # np.testing.assert_array_equal(f.update([1.], [1000.]), [10.0])

        f = pg.trans.TransCumulative()
        f.add(pg.trans.TransLog(), 5)
        f.add(pg.trans.TransLog(), 5)

        np.testing.assert_array_equal(f.at(0).fwd(np.ones(10)*10),
                                      np.log(np.ones(10)*10))
        np.testing.assert_array_equal(f.fwd(np.ones(10)*10),
                                      np.log(np.ones(10)*10))
        # tm2 = pg.trans.TransLog()
        # tc.add(tm2, 5, 10)
        # fop._modelTrans = pg.trans.TransCumulative()
        # fop._modelTrans.add(tm2, size=nModel)


        #fop._modelTrans = pg.trans.TransLog()


    def test_DataContainerFilter(self):
        """
        """
        data = pg.DataContainer()
        data.resize(5)

        data.markValid([0, 4])
        self.assertEqual(data('valid'), [1.0, 0.0, 0.0, 0.0, 1.0])

        data.markInvalid(pg.core.IndexArray(np.arange(5, dtype="long")))
        self.assertEqual(data('valid'), [0.0, 0.0, 0.0, 0.0, 0.0])

        data.markValid(np.arange(5, dtype="long"))
        self.assertEqual(data('valid'), [1.0, 1.0, 1.0, 1.0, 1.0])

        data.markInvalid(range(5))
        self.assertEqual(data('valid'), [0.0, 0.0, 0.0, 0.0, 0.0])

        x = np.arange(5, dtype='float')

        data.markValid(pg.Vector(x) > 2.0)
        self.assertEqual(data('valid'), [0.0, 0.0, 0.0, 1.0, 1.0])

        data.markValid(pg.BVector(x < 2.0))
        self.assertEqual(data('valid'), [1.0, 1.0, 0.0, 1.0, 1.0])

        data.markInvalid(pg.find(x > 3.0))
        self.assertEqual(data('valid'), [1.0, 1.0, 0.0, 1.0, 0.0])

        data.markInvalid(x < 1.0)
        self.assertEqual(data('valid'), [0.0, 1.0, 0.0, 1.0, 0.0])


    def test_DataContainerSensors(self):
        data = pg.DataContainer()

        sensors = [[x, 0.0] for x in range(5)]
        data.setSensorPositions(sensors)
        data.setSensorPositions(data.sensors()[::-1])

        self.assertEqual(data.sensor(0), [4., 0.0, 0.0])
        self.assertEqual(data.sensor(4), [0., 0.0, 0.0])


    def test_DataContainerIndex(self):
        data = pg.DataContainer()
        data['b'] = np.ones(2) * 3.14
        np.testing.assert_array_equal(data['b'], np.ones(2)*3.14)
        self.assertEqual(type(data['b']), type(pg.Vector()))

        data['b'][0] = 1.0
        self.assertEqual(data['b'][0], 1.0)

        data.registerSensorIndex('a')
        data['a'] = np.ones(2)
        np.testing.assert_array_equal(data['a'], np.ones(2))
        self.assertEqual(type(data['a']), type(np.array(1)))
        self.assertEqual(data['a'].dtype, 'int')
        data['a'][0] = 1.0 # will not work for sensorIndex until its changed in the datacontainer as IndexArray

        data['a'] = np.ones(2)*1.2
        np.testing.assert_array_equal(data['a'], np.ones(2))
        self.assertEqual(type(data['a']), type(np.array(1)))
        self.assertEqual(data['a'].dtype, 'int')


    def test_Operators(self):
        t = pg.Vector(10, 1.0)
        self.assertEqual(len(t == 1.0), len(t > 0))
        self.assertEqual(len(t == 1.0), len(t == 1))


    def test_Int64Problem(self):
        data = pg.DataContainerERT()
        data.createFourPointData(0, 0, 1, 2, 3)
        pos = np.arange(4, dtype=int)
        data.createFourPointData(1, pos[0], pos[1], pos[2], pos[3])
        pos = np.arange(4, dtype=np.int32)
        data.createFourPointData(2, pos[0], pos[1], pos[2], pos[3])
        pos = np.arange(4, dtype=np.int64)
        data.createFourPointData(3, pos[0], pos[1], pos[2], pos[3])
        pos = np.arange(4, dtype=float)
        data.createFourPointData(4, pos[0], pos[1], pos[2], pos[3])
        pos = np.arange(4, dtype=np.float32)
        data.createFourPointData(5, pos[0], pos[1], pos[2], pos[3])
        pos = np.arange(4, dtype=np.float64)
        data.createFourPointData(6, pos[0], pos[1], pos[2], pos[3])
        pos = np.arange(4)
        data.createFourPointData(7, pos[0], pos[1], pos[2], pos[3])
        pos = range(4)
        data.addFourPointData(pos[0], pos[1], pos[2], pos[3])
        #print(data('a'), data('b'), data('m'), data('n'))
        self.assertEqual(sum(data('a')), 9*0)
        self.assertEqual(sum(data('b')), 9*1)
        self.assertEqual(sum(data('m')), 9*2)
        self.assertEqual(sum(data('n')), 9*3)


    def test_PosConstMember(self):
        p1 = pg.Pos(1.0, 0.0, 0.0)
        p2 = pg.Pos(0.0, 1.0, 0.0)

        p3 = p1.cross(p2)
        self.assertEqual(p3, pg.Pos(0.0, 0.0, 1.0))


    def test_Hash(self):
        v1 = pg.Vector(10, 2.)
        v2 = pg.Vector(10, 2.)

        self.assertFalse(pg.Vector(1, 0.).hash() == pg.Vector(2, 0.).hash())

        self.assertEqual(v1.hash(), v2.hash())
        self.assertEqual(hash(v1), hash(v2))
        v2[2] = 3.
        self.assertFalse(v1.hash() is v2.hash())
        v2[2] = 2.
        self.assertTrue(v1.hash() == v2.hash())
        self.assertEqual(v1.hash(), pg.Vector(10, 2.).hash())

        self.assertEqual(pg.utils.valHash([1, 0, 0])==
                         pg.utils.valHash([0, 1, 0]), False)


    def test_HashData(self):
        d1 = pg.DataContainerERT()
        d2 = pg.DataContainerERT()

        self.assertEqual(d1.hash(), d2.hash())
        d1.createSensor([1.0, 0.0])
        d2.createSensor([2.0, 0.0])
        self.assertFalse(d1.hash() == d2.hash())
        d2.setSensor(0, [1.0, 0.0])
        self.assertTrue(d1.hash() == d2.hash())

        d1.resize(10)
        d2.resize(12)
        d1.add('a', pg.Vector(d1.size(), 1.0))
        d2.add('a', pg.Vector(d2.size(), 1.0))
        self.assertFalse(d1.hash() == d2.hash())

        d2.resize(10)
        self.assertTrue(d1.hash() == d2.hash())
        d2('a')[3] = 2.0
        self.assertFalse(d1.hash() == d2.hash())
        d2('a')[3] = 1.0
        self.assertTrue(d1.hash() == d2.hash())

    def test_HashMesh(self):
        m1 = pg.Mesh()
        m2 = pg.Mesh()

        self.assertTrue(m1.hash() == m2.hash())

        m1.createNode([1.0, 0.0])
        m2.createNode([2.0, 0.0])
        self.assertFalse(m1.hash() == m2.hash())
        m2.node(0).setPos([1.0, 0.0])
        self.assertTrue(m1.hash() == m2.hash())


    # does not work .. need time to implement
    # def test_DataContainerWrite(self):
    #     data = pg.DataContainer()
    #     data.save('test.dat')
    #     fi = open('test2.dat', 'w')
    #     data.write(fi)
    #     fi.close()

    def test_DataTypes(self):
        pg.core.showSizes()


if __name__ == '__main__':
    pg.core.setDeepDebug(0)
    unittest.main()
