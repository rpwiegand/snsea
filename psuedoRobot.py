import numpy as np

SHOWPLOTS = True

if SHOWPLOTS:
  import matplotlib.pyplot as plt


class BaseController:
    def getBehaviorVector(self, controlVector):
        return np.array(controlVector)


class PsuedoRobotController(BaseController):
    def __init__(self, maxAngleRadians, maxLength):
        self.maxAngleRadians = maxAngleRadians
        self.maxLength = maxLength
        self.wayPoints = []

    def getBehaviorVector(self, controlVector):
        numPoints = len(controlVector)
        perPointDist = self.maxLength / numPoints


        # Make sure all values are between -1 and 1:
        #maxVal = np.max(np.abs(controlVector))
        #if maxVal > 1:
        #    controlVector /= maxVal

        # Set initial orientation based on the first control angle between -pi and pi
        orientation = controlVector[0] * np.pi

        # Set the rest of the angles between +/- maxAngle
        controlVector = controlVector[1:-1] * 1 #self.maxAngleRadians  ## RPW:  Maybe don't do this ... let the turn angle be as severe as you like

        self.wayPoints = [ (0,0) ]
        for controlAngle in controlVector:
            x = self.wayPoints[-1][0]
            y = self.wayPoints[-1][1]

            orientation += controlAngle
            x += perPointDist * np.cos(orientation)
            y += perPointDist * np.sin(orientation)

            self.wayPoints.append( (x,y) )

        return np.array(self.wayPoints[-1])


    def printWayPoints(self):
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



