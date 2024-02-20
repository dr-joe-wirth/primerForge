import time, unittest
from bin.Clock import Clock

class ClockTest(unittest.TestCase):
    """evaluating the Clock class
    """
    # constant
    SLEEP_TIME = 4
    
    def setUp(self) -> None:
        self.clock = Clock()
        self.start = time.time()
    
    def tearDown(self) -> None:
        return super().tearDown()
    
    def _timeToSeconds(hrs:int, min:int, sec:float) -> float:
        """converts hrs, min, and sec to a seconds

        Args:
            hrs (int): hours
            min (int): minutes
            sec (float): seconds

        Returns:
            float: total time in seconds
        """
        return hrs*3600 + min*60 + sec
    
    def testA_getTimeFormat(self) -> None:
        """is Clock.getTime returning the expected format
        """
        clockTime = self.clock.getTime()
        
        self.assertIsInstance(clockTime, tuple)
        self.assertEqual(len(clockTime), 3)
        self.assertIsInstance(clockTime[0], int)
        self.assertIsInstance(clockTime[1], int)
        self.assertIsInstance(clockTime[2], float)
    
    def testB_restart(self) -> None:
        """does Clock.restart work
        """
        # restart the clock then get the time
        self.clock.restart()
        hrs,min,sec = self.clock.getTime()
        
        # make sure almost 0 time has elapsed
        self.assertEqual(hrs, 0)
        self.assertEqual(min, 0)
        self.assertEqual(round(sec), 0)
        
        # restart the clock, sleep, then get the time again
        self.clock.restart()
        time.sleep(ClockTest.SLEEP_TIME)
        hrs,min,sec = self.clock.getTime()
        
        # make sure the time matches expected value
        self.assertEqual(hrs, 0)
        self.assertEqual(min, 0)
        self.assertEqual(round(sec), ClockTest.SLEEP_TIME)
    
    def testC_getTimeAccuracy(self) -> None:
        """is Clock.getTime accurate
        """
        # sleep for a few seconds
        time.sleep(ClockTest.SLEEP_TIME)
        
        # get the clock time the test duration
        hrs,min,sec = self.clock.getTime()
        testTime = time.time() - self.start

        # convert the clock time to seconds
        clockTime = ClockTest._timeToSeconds(hrs, min, sec)
    
        # make sure the durations are nearly identical
        self.assertAlmostEqual(clockTime, testTime, places=1)
    
    def testD_getTimeString(self):
        """is Clock.getTimeString accurate
        """
        # constant
        NUM_DEC = 3
        
        # sleep for a few seconds
        time.sleep(ClockTest.SLEEP_TIME)
        
        # extract times and extract the string
        testTime = time.time() - self.start
        hrs,min,sec = self.clock.getTime(NUM_DEC)
        timeStr = self.clock.getTimeString(NUM_DEC)
        
        # make sure it is a string
        self.assertIsInstance(timeStr, str)
        
        # extract individual values from the string
        h,m,s = timeStr.split(":")
        h = int(h)
        m = int(m)
        s = float(s)
        
        # make sure the values look appropriate
        self.assertEqual(hrs, h)
        self.assertEqual(min, m)
        self.assertAlmostEqual(sec, s, places=NUM_DEC)
        
        # make sure the time string looks right
        self.assertEqual(timeStr, f"00:00:0{sec:.{NUM_DEC}f}")
        self.assertEqual(timeStr, f"00:00:0{testTime:.{NUM_DEC}f}")
