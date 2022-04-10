#include <iostream>
#include <math.h>
#include <string.h>

#ifndef DEBUG
#define DEBUG false
#endif

using namespace std;

class GradientDescent {
public:
    int nDim;
    double* w;
    double alpha;
    double tolerance;
    bool isMin;
    int nLimits;
    double** limits;

    GradientDescent(int nDim, double* w, double alpha, double tolerance, bool isMin) {
        this->nDim = nDim;
        this->w = w;
        this->alpha = alpha;
        this->tolerance = tolerance;
        this->isMin = isMin;
    }

    GradientDescent(int nDim, double* w) {
        this->nDim = nDim;
        this->w = w;
        this->alpha = 0.1;
        this->tolerance = 0.01;
        this->isMin = true;
    }

    ~GradientDescent() {
        delete[] this->w;
        delete[] this->limits;
    }

    void printPoint(double* point, int nDim) {
        int i = 0;
        cout << "(";
        for (i; i < nDim; ++i) {
            cout << point[i];
            cout << ((i < nDim - 1) ? ", " : ")");
        }
        cout << endl;
    }

    double* getGradient() {
        double* gradient = new double[this->nDim - 1];
        copy(this->w, this->w + this->nDim - 1, gradient);
        return gradient;
    }

    void setBorders(int n, double** limits) {
        this->nLimits = n;
        this->limits = limits;
    }

    bool isInLimits(double* x) {
        for (int i = 0; i < this->nLimits; ++i) {
            double result = 0;
            for (int j = 0; j < this->nDim; ++j) {
                if (j < this->nDim - 1) {
                    result += this->limits[i][j] * floor(x[j] * 10000000.0) / 10000000.0; // xi case except y (x2)
                } else {
                    result += this->limits[i][j]; // independent case
                }
            }
            if (result < 0) {
                return false;
            }
        }
        return true;
    }

    void getPerpendicular(double* point, double* direction, double* perpendicular) {
        copy(direction, direction + this->nDim - 1, perpendicular);
        perpendicular[0] = -perpendicular[0] / perpendicular[1];
        perpendicular[1] = -1;
        perpendicular[2] = perpendicular[0] * (-point[0]) + point[1];
    }

    bool getIntersectionPoint(double* line1, double* line2, double*& intersectionPoint) {
        double* temp = new double[this->nDim];
        double factor = 0;
        if (line1[0] == 0 && line2[0] == 0) {
            delete[] temp;
            delete[] intersectionPoint;
            intersectionPoint = nullptr;
            if (line1[1] == 0 && line2[1] == 0) { // all = 0
                cout << "ERROR: Invalid lines in getIntersectionPoint function" << endl;
                return false;
            }
            if (line1[1] == line2[1] && line1[2] == line2[2]) { // the same lines
                return true;
            }
            // parallel lines
            return false;
        }
        if (line1[0] == 0) {
            for (int i = 0; i < this->nDim; ++i) {
                temp[i] = line1[i];
            }
        } else if (line2[0] == 0) {
            for (int i = 0; i < this->nDim; ++i) {
                temp[i] = line2[i];
            }
        } else {
            factor = -(line1[0] / line2[0]);
            for (int i = 0; i < this->nDim; ++i) {
                temp[i] = (factor * line2[i]) + line1[i];
            }
        }
        if (temp[1] != 0 && temp[1] != INFINITY) {
            intersectionPoint[1] = -temp[2] / temp[1];
            intersectionPoint[0] = -(line1[1] * intersectionPoint[1] + line1[2]) / line1[0];
            delete[] temp;
            return true;
        } else {
            bool isIntersect = temp[2] == 0 || isnan(temp[2]);
            delete[] temp;
            delete[] intersectionPoint;
            intersectionPoint = nullptr;
            return isIntersect;
        }
    }

    void setMidpoint(double* point, double* firstIntersectionPoint, double* secondIntersectionPoint) {
        for (int i = 0; i < nDim - 1; ++i) {
            point[i] = (firstIntersectionPoint[i] + secondIntersectionPoint[i]) / 2;
        }
#if (DEBUG)
        cout << "First intersection point: ";
        printPoint(firstIntersectionPoint, this->nDim - 1);
        cout << "Second intersection point: ";
        printPoint(secondIntersectionPoint, this->nDim - 1);
        cout << "New point (midpoint): ";
        printPoint(point, this->nDim - 1);
#endif
    }

    double distance(double* point1, double* point2, int nDim) {
        double distance = 0;
        for (int i = 0; i < nDim - 1; ++i) {
            distance += pow(point1[i] - point2[i], 2);
        }
        return sqrt(distance);
    }

    void move(double* point, double* direction) {
        for (int i = 0; i < nDim - 1; ++i) {
            point[i] = point[i] + this->alpha * direction[i];
        }
    }

    double* run(int maxSteps) {
        int step = 1;
        int i = 0;
        double* point = new double[this->nDim - 1];
        double* prevPoint = new double[this->nDim - 1];
        double* gradient = getGradient();
        double* perpendicular = new double[this->nDim - 1];
        bool isIntersect = false;
        double* firstIntersectionPoint = nullptr;
        double* secondIntersectionPoint = nullptr;
        bool extremumPointFound = false;
        bool isPrevPointExtremum = false;

        if (this->isMin) {
            for (i = 0; i < nDim - 1; ++i) {
                gradient[i] = -gradient[i];
            }
        }

        for (i = 0; i < nDim - 1; ++i) {
            point[i] = 0;
            prevPoint[i] = 0;
        }
#if (DEBUG)
        cout << "Starting point: ";
        printPoint(point, this->nDim - 1);
#endif

        while (step <= maxSteps && !extremumPointFound) {
#if (DEBUG)
            cout << "---------" << endl;
            cout << "Step N" << step << endl;
#endif

            for (i = 0; i < nDim - 1; ++i) {
                prevPoint[i] = point[i];
            }
            move(point, gradient);
#if (DEBUG)
            cout << "x = ";
            printPoint(point, this->nDim - 1);
#endif

            if (distance(point, prevPoint, this->nDim - 1) < this->tolerance) {
                cout << "Stopped due to tolerance." << endl;
                extremumPointFound = true;
                break;
            }
            if (!isInLimits(point)) {
                getPerpendicular(prevPoint, gradient, perpendicular);
#if (DEBUG)
                cout << "Out of borders" << endl;
                cout << "perpendicular = ";
                printPoint(perpendicular, this->nDim);
#endif
                for (i = 0; i < this->nLimits; ++i) {
                    double* intersectionPoint = new double[this->nDim - 1];
                    isIntersect = getIntersectionPoint(perpendicular, this->limits[i], intersectionPoint);
                    if (isIntersect) {
                        if (intersectionPoint == nullptr) {
                            cout << this->limits[i][0] << this->limits[i][1] << this->limits[i][2] << endl;
                            cout << "The extremum points are on the N" << i + 1 << " border." << endl;
                            extremumPointFound = true;
                            isPrevPointExtremum = true;
                            break;
                        }
                        if (!isInLimits(intersectionPoint)) {
                            isIntersect = false;
                            delete[] intersectionPoint;
                            intersectionPoint = nullptr;
                        } else if (firstIntersectionPoint == nullptr) {
                            firstIntersectionPoint = intersectionPoint;
                        } else {
                            secondIntersectionPoint = intersectionPoint;
                            setMidpoint(point, firstIntersectionPoint, secondIntersectionPoint);
                            break;
                        }
                    }
                }
                if (firstIntersectionPoint != nullptr) {
                    if (secondIntersectionPoint != nullptr && distance(firstIntersectionPoint, secondIntersectionPoint, this->nDim - 1) > this->tolerance) {
                        delete[] secondIntersectionPoint;
                        secondIntersectionPoint = nullptr;
                    } else { // Only 1 intersection point was found
                        extremumPointFound = true;
                        isPrevPointExtremum = true;
                        delete[] firstIntersectionPoint;
                        firstIntersectionPoint = nullptr;
                        break;
                    }
                    delete[] firstIntersectionPoint;
                    firstIntersectionPoint = nullptr;
                } else {
                    cout << "No intersection was found. It means that the previous point was out of limits too." << endl;
                    extremumPointFound = true;
                    isPrevPointExtremum = true;
                    break;
                }

                if (distance(point, prevPoint, this->nDim - 1) < this->tolerance) {
#if (DEBUG)
                    cout << "New point is too close to the mid point. Decreasing the alpha." << endl;
#endif
                    this->alpha /= 2;
                }
            }
            step++;
        }
        if (isPrevPointExtremum) {
            for (i = 0; i < nDim - 1; ++i) {
                point[i] = prevPoint[i];
            }
        }
        cout << "Steps: " << step - 1 << " (Max steps: " << maxSteps << ")" << endl;
        delete[] prevPoint;
        delete[] gradient;
        return point;
    }
};

int main() {
    int nDim;
    double* w;
    int nLimits;
    double** limits;
    double alpha;
    double tolerance;
    char isMinStr[3];
    bool isMin;
    int maxSteps;
    int i = 0;
    int j = 0;
    double* point;

    cout << "Enter the number of dimensions: ";
    cin >> nDim;
    w = new double[nDim];
    for (i = 0; i < nDim; ++i) {
        if (i < nDim - 1) {
            cout << "w[" << j << "] = ";
        } else {
            cout << "b = ";
        }
        cin >> w[i];
    }

    cout << "Enter the count of limits: ";
    cin >> nLimits;
    limits = new double*[nLimits];
    for (i = 0; i < nLimits; ++i) {
        limits[i] = new double[nDim];
        cout << "Enter limit N" << i << ":" << endl;
        for (j = 0; j < nDim; ++j) {
            if (j < nDim - 1) {
                cout << "w[" << j << "] = ";
            } else {
                cout << "b = ";
            }
            cin >> limits[i][j];
        }
        limits[i];
    }

    cout << "Enter min or max: ";
    cin >> isMinStr;
    isMin = strcmp(isMinStr, "max") == 0 ? false : true;
    cout << "Enter the alpha: ";
    cin >> alpha;
    cout << "Enter the tolerance: ";
    cin >> tolerance;
    cout << "Enter the maximum of steps count: ";
    cin >> maxSteps;

    GradientDescent* gradientDescent = new GradientDescent(nDim, w, alpha, tolerance, isMin);
    gradientDescent->setBorders(nLimits, limits);
    cout << "Running..." << endl;
    point = gradientDescent->run(maxSteps);
    cout << "Extremum point: ";
    gradientDescent->printPoint(point, nDim - 1);

    delete[] point;
    delete[] gradientDescent;
    cout << "End" << endl;
    return 0;
}
