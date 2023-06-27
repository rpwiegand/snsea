## This is a very simple robot controller that takes a control vector and returns a path
## of way points based on that control vector.  This is intended to experiment with the 
## mapping of genotype to behavior in terms of novelty search.
import numpy as np

SHOWPLOTS = True  # RPW:  If true, the unit test will use Matplotlib to show the behavior

if SHOWPLOTS:
  import matplotlib.pyplot as plt


class BaseController:
    def getBehaviorVector(self, controlVector):
        return np.array(controlVector)


class PsuedoRobotController(BaseController):
    def __init__(self, maxAngleRadians, maxLength):
        """
        Contructor for the robot controller.  This sets up the maximum turn agnel the the
        maximum length of the total path.
        """
        self.maxAngleRadians = maxAngleRadians
        self.maxLength = maxLength
        self.wayPoints = []


    def getBehaviorVector(self, controlVector):
        """
        Given the control vector, return the way points of the robot.  The first number in the control
        vector is the initial orientation.  The remaining numbers are angles (in radians) of the turn.
        """
        # Estimate the distance between each point
        numPoints = len(controlVector)
        perPointDist = self.maxLength / numPoints


        # Set initial orientation based on the first control angle between -pi and pi
        orientation = controlVector[0] * np.pi

        # Set the rest of the angles between +/- maxAngle
        controlVector = controlVector[1:-1] * 1 #RPW:  Changed my mind ... let the angles range as much as they want ...

        self.wayPoints = [ (0,0) ]
        for controlAngle in controlVector:
            x = self.wayPoints[-1][0]
            y = self.wayPoints[-1][1]

            # Turn, take a step, and update the position
            orientation += controlAngle
            x += perPointDist * np.cos(orientation)
            y += perPointDist * np.sin(orientation)

            # Store the position in the waypoints vector
            self.wayPoints.append( (x,y) )

        return np.array(self.wayPoints[-1])


    def printWayPoints(self):
        """
        Print the robot's way points
        """
        print()
        print("X", "\t", "Y")
        for pt in self.wayPoints:
            print(pt[0], "\t", pt[1])
        print()
    


def unitTest():
    controller = PsuedoRobotController(np.pi/3,10)

    for trial in range(10):
      randomTurnings = np.random.uniform(low=-1.0, high=1.0, size=11)
      controller.getBehaviorVector(randomTurnings)
      controller.printWayPoints()

      if SHOWPLOTS:
        x, y = zip(*controller.wayPoints)
        plt.plot(x,y)

    if SHOWPLOTS:
      plt.xlim([-10,10])
      plt.ylim([-10,10])
      plt.show()


    

if __name__ == "__main__":
    unitTest()



