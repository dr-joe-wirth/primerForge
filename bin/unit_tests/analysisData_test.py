from __future__ import annotations
import unittest
from bin.Primer import Primer
from itertools import combinations
from bin.AnalysisData import AnalysisData, Level


class AnalysisDataTest(unittest.TestCase):
    # global attributes
    LEVEL_STR = ('candidate kmer', 'unfiltered pair', 'filtered pair', 'final pair')
    PRIMER = Primer('atcg', 'contig', 8, 4, Primer.PLUS)
    INDEX = 0
    NAME = 'genome' 
    
    def setUp(self) -> None:
        return super().setUp()
    
    def tearDown(self) -> None:
        return super().tearDown()
    
    def testA_constructor(self) -> None:
        """make sure attributes have expected values
        """
        data = AnalysisData(AnalysisDataTest.PRIMER, AnalysisDataTest.INDEX, AnalysisDataTest.NAME)
        self.assertEqual(data.primer, self.PRIMER)
        self.assertEqual(data.getIndex(), self.INDEX)
        self.assertEqual(data.name, self.NAME)
        
    def testB_setGetLevel(self) -> None:
        """checks that Level.setLevel, AnalysisData.setLevel, and AnalysisData.getLevel are working
        """
        # initialize a Level object
        level = Level()
        
        # make sure the constructor worked
        self.assertEqual(str(level), AnalysisDataTest.LEVEL_STR[0])
        
        # check that setLevel and str works for all levels
        for lev in AnalysisDataTest.LEVEL_STR:
            # make sure setLevel works with a string
            level.setLevel(lev)
            self.assertEqual(str(level), lev)
            
            # make sure setLevel works with a Level
            newLev = Level()
            newLev.setLevel(lev)
            level.setLevel(newLev)
            self.assertEqual(str(level), lev)
        
        # check that setLevel fails with invalid inputs
        self.assertRaises(ValueError, level.setLevel, 'bad level')
        self.assertRaises(ValueError, level.setLevel, 0)

        # initialize an AnalysisData object
        data = AnalysisData(AnalysisDataTest.PRIMER, AnalysisDataTest.INDEX, AnalysisDataTest.NAME)
        
        # check that setLevel and getLevel works for all levels
        for lev in AnalysisDataTest.LEVEL_STR:
            # make sure setLevel works with a string
            data.setLevel(lev)
            self.assertEqual(data.getLevel(), lev)
            
            # make sure setLevel works with a Level
            newLev = Level()
            newLev.setLevel(lev)
            level.setLevel(newLev)
            self.assertEqual(data.getLevel(), lev)
        
        # check that setLevel fails with invalid inputs
        self.assertRaises(ValueError, data.setLevel, 'bad level')
        self.assertRaises(ValueError, data.setLevel, 0)
    
    def testC_LevelEquality(self) -> None:
        """checks that the == operator works for Level class
        """
        # initialize levels
        level = Level()
        other = Level()
        
        # they should be equal
        self.assertEqual(level, other)
        
        # make sure they stay equal for each level
        for lev in AnalysisDataTest.LEVEL_STR:
            level.setLevel(lev)
            other.setLevel(lev)
            
            # should work for levels
            self.assertEqual(level, other)
            
            # should work for strings
            self.assertEqual(level, lev)
        
        # make sure it fails with int
        self.assertRaises(TypeError, level.__eq__, 1)
        
    def testD_LevelInequality(self) -> None:
        """checks that != operator works for Level class
        """
        # initialize variables
        level = Level()
        other = Level()
        
        # make sure all combinations of levels are not equal
        for lev1,lev2 in combinations(AnalysisDataTest.LEVEL_STR, 2):
            level.setLevel(lev1)
            other.setLevel(lev2)
            
            # should work for levels
            self.assertNotEqual(level, other)
            
            # should work for strings
            self.assertNotEqual(level, lev2)
            self.assertNotEqual(other, lev1)
        
        # make sure it fails with int
        self.assertRaises(TypeError, level.__ne__, 1)
    
    def testE_LevelAdd(self) -> None:
        """makes sure the + operator works for Level class
        """
        # for each level
        for idx in range(len(AnalysisDataTest.LEVEL_STR)):
            # add to the level
            level = Level() + idx
            
            # create a level at the same level
            other = Level()
            other.setLevel(AnalysisDataTest.LEVEL_STR[idx])
            
            # make sure they are equal
            self.assertEqual(level, other)
        
        # make sure cannot add beyond bounds
        level = Level()
        self.assertRaises(Exception, level.__add__, -1)
        self.assertRaises(Exception, level.__add__, 4)
        
        # make sure cannot add non ints
        level = Level()
        self.assertRaises(TypeError, level.__add__, 1.2)
        self.assertRaises(TypeError, level.__add__, Level())
    
    def testF_LevelSelfAdd(self) -> None:
        """checks that the += operator works for Level class
        """
        # for each level
        for idx in range(len(AnalysisDataTest.LEVEL_STR)):
            # initialize levels at the floor
            level = Level()
            other = Level()
            
            # set the levels both ways
            level += idx
            other.setLevel(AnalysisDataTest.LEVEL_STR[idx])
            
            # make sure they are equal
            self.assertEqual(level, other)
        
        # make sure cannot add beyond bounds
        level = Level()
        self.assertRaises(Exception, level.__iadd__, -1)
        self.assertRaises(Exception, level.__iadd__, 4)
        
        # make sure cannot add non ints
        level = Level()
        self.assertRaises(TypeError, level.__iadd__, 1.2)
        self.assertRaises(TypeError, level.__iadd__, Level())

    def testG_LevelSub(self) -> None:
        """checks that - operator works for Level class
        """
        # subtract numbers to bring a higher level back to floor
        for idx in range(len(AnalysisDataTest.LEVEL_STR)-1,-1,-1):
            # create level at ceiling, then subtract to the floor
            level = Level()
            level.setLevel(AnalysisDataTest.LEVEL_STR[idx])
            level = level - idx
            
            # should be equal
            self.assertEqual(level, Level())
        
        # subtract 1 
        for idx in range(1,len(AnalysisDataTest.LEVEL_STR)):
            # initialize variables
            level = Level()
            other = Level()
            
            # make the level one less
            level.setLevel(AnalysisDataTest.LEVEL_STR[idx])
            level = level - 1
            
            # set the other level to one less
            other.setLevel(AnalysisDataTest.LEVEL_STR[idx-1])
            
            # they should be equal
            self.assertEqual(level, other)
        
        # make sure cannot sub beyond bounds
        level = Level()
        self.assertRaises(Exception, level.__sub__, 1)
        self.assertRaises(Exception, level.__sub__, -4)
        
        # make sure cannot sub non ints
        level = Level()
        level.setLevel(AnalysisDataTest.LEVEL_STR[-1])
        self.assertRaises(TypeError, level.__sub__, 1.2)
        self.assertRaises(TypeError, level.__sub__, Level())
    
    def testH_LevelSelfSub(self) -> None:
        """checks that -= operator works for Level class
        """
        # subtract numbers to bring a higher level back to floor
        for idx in range(len(AnalysisDataTest.LEVEL_STR)-1,-1,-1):
            # create level at ceiling, then subtract to the floor
            level = Level()
            level.setLevel(AnalysisDataTest.LEVEL_STR[idx])
            level -= idx
            
            # should be equal
            self.assertEqual(level, Level())
        
        # subtract 1 
        for idx in range(1,len(AnalysisDataTest.LEVEL_STR)):
            # initialize variables
            level = Level()
            other = Level()
            
            # make the level one less
            level.setLevel(AnalysisDataTest.LEVEL_STR[idx])
            level -= 1
            
            # set the other level to one less
            other.setLevel(AnalysisDataTest.LEVEL_STR[idx-1])
            
            # they should be equal
            self.assertEqual(level, other)
        
        # make sure cannot sub beyond bounds
        level = Level()
        self.assertRaises(Exception, level.__isub__, 1)
        self.assertRaises(Exception, level.__isub__, -4)
        
        # make sure cannot sub non ints
        level = Level()
        level.setLevel(AnalysisDataTest.LEVEL_STR[-1])
        self.assertRaises(TypeError, level.__isub__, 1.2)
        self.assertRaises(TypeError, level.__isub__, Level())
    
    def testI_LevelGreaterThan(self) -> None:
        """check that the > operator works for Level class
        """
        # for each level
        for idx in range(1,len(AnalysisDataTest.LEVEL_STR)):
            # create a new level at that level
            level = Level() + idx
            
            # make sure it is greater than the floor and all others below it
            for less in range(idx):
                other = Level() + less
                
                # should work with levels
                self.assertGreater(level, other)
                
                # should work with strings
                self.assertGreater(level, AnalysisDataTest.LEVEL_STR[less])
        
        # make sure it fails on int
        self.assertRaises(TypeError, level.__gt__, 1)
    
    def testJ_LevelLessThan(self) -> None:
        """check that < operator works for Level class
        """
        # for each level
        for idx in range(1,len(AnalysisDataTest.LEVEL_STR)):
            # create a new level at that level
            level = Level() + idx
            
            # make sure it is less than all other above it
            for less in range(idx):
                other = Level() + less
                self.assertLess(other, level)
            
            # make sure the floor is always less
            self.assertLess(Level(), level)

        # make sure it fails on int
        self.assertRaises(TypeError, level.__lt__, 1)
        
    def testK_LevelGreaterThanEqual(self) -> None:
        """checks that >= operator works for Level class
        """
        # for each level
        for idx in range(len(AnalysisDataTest.LEVEL_STR)):
            # create a new level at that level
            level = Level() + idx
            
            # make sure it is greater than or equal to all others below it
            for less in range(idx):
                other = Level() + less
                self.assertGreaterEqual(level, other)
                
                other.setLevel(str(level))
                self.assertGreaterEqual(level, other)
                
            # make sure it is greater than or equal to the floor
            self.assertGreaterEqual(level, Level())
        
        # make sure it fails on int
        self.assertRaises(TypeError, level.__ge__, 1)
    
    def testL_LevelLessThanEqual(self) -> None:
        """checks that <= operator works for Level class
        """
        # for each level
        for idx in range(len(AnalysisDataTest.LEVEL_STR)):
            # create a new level at that level
            level = Level() + idx
            
            # make sure it is less than or equal to all others above it
            for less in range(idx):
                other = Level() + less
                self.assertLessEqual(other, level)
                
                other.setLevel(str(level))
                self.assertLessEqual(other, level)
                
            # make sure it is greater than or equal to the floor
            self.assertLessEqual(Level(), level)
        
        # make sure it fails on int
        self.assertRaises(TypeError, level.__le__, 1)
    
    def testM_incrementLevel(self) -> None:
        """check that AnalysisData.incrementLevel works
        """
        # create objects
        data = AnalysisData(AnalysisDataTest.PRIMER, AnalysisDataTest.INDEX, AnalysisDataTest.NAME)
        level = Level()
        
        # for each level above the floor
        for idx in range(1,len(AnalysisDataTest.LEVEL_STR)):
            # increase the level and increment data's level
            level += 1
            data.incrementLevel()
            
            # they should be equal
            self.assertEqual(data.getLevel(), level)
        
        # make sure it fails when at the ceiling
        self.assertRaises(Exception, data.incrementLevel)
    
    def testN_getUpdatePairs(self) -> None:
        """checks that AnalysisData.getPairs and AnalysisData.updatePairs work
        """
        # create object
        data = AnalysisData(AnalysisDataTest.PRIMER, AnalysisDataTest.INDEX, AnalysisDataTest.NAME)
        
        # should be an empty list
        self.assertIsInstance(data.getPairs(), list)
        self.assertListEqual(data.getPairs(), list())
        
        # updating should add to the list
        data.updatePairs(1)
        self.assertListEqual(data.getPairs(), [1])
        
        # updating should add to the list
        data.updatePairs(2)
        self.assertListEqual(data.getPairs(), [1, 2])
        
        # redundant additions should not be present
        data.updatePairs(2)
        self.assertListEqual(data.getPairs(), [1, 2])
