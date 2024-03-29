//=============================================================================================
// Szamitogepes grafika hazi feladat keret. Ervenyes 2014-tol.          
// A //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// sorokon beluli reszben celszeru garazdalkodni, mert a tobbit ugyis toroljuk. 
// A beadott program csak ebben a fajlban lehet, a fajl 1 byte-os ASCII karaktereket tartalmazhat. 
// Tilos:
// - mast "beincludolni", illetve mas konyvtarat hasznalni
// - faljmuveleteket vegezni (printf is fajlmuvelet!)
// - new operatort hivni az onInitialization függvényt kivéve, a lefoglalt adat korrekt felszabadítása nélkül 
// - felesleges programsorokat a beadott programban hagyni
// - tovabbi kommenteket a beadott programba irni a forrasmegjelolest kommentjeit kiveve
// ---------------------------------------------------------------------------------------------
// A feladatot ANSI C++ nyelvu forditoprogrammal ellenorizzuk, a Visual Studio-hoz kepesti elteresekrol
// es a leggyakoribb hibakrol (pl. ideiglenes objektumot nem lehet referencia tipusnak ertekul adni)
// a hazibeado portal ad egy osszefoglalot.
// ---------------------------------------------------------------------------------------------
// A feladatmegoldasokban csak olyan gl/glu/glut fuggvenyek hasznalhatok, amelyek
// 1. Az oran a feladatkiadasig elhangzottak ES (logikai AND muvelet)
// 2. Az alabbi listaban szerepelnek:  
// Rendering pass: glBegin, glVertex[2|3]f, glColor3f, glNormal3f, glTexCoord2f, glEnd, glDrawPixels
// Transzformaciok: glViewport, glMatrixMode, glLoadIdentity, glMultMatrixf, gluOrtho2D, 
// glTranslatef, glRotatef, glScalef, gluLookAt, gluPerspective, glPushMatrix, glPopMatrix,
// Illuminacio: glMaterialfv, glMaterialfv, glMaterialf, glLightfv
// Texturazas: glGenTextures, glBindTexture, glTexParameteri, glTexImage2D, glTexEnvi, 
// Pipeline vezerles: glShadeModel, glEnable/Disable a kovetkezokre:
// GL_LIGHTING, GL_NORMALIZE, GL_DEPTH_TEST, GL_CULL_FACE, GL_TEXTURE_2D, GL_BLEND, GL_LIGHT[0..7]
//
// NYILATKOZAT
// ---------------------------------------------------------------------------------------------
// Nev    : Boczán Tamás
// Neptun : A5X61F
// ---------------------------------------------------------------------------------------------
// ezennel kijelentem, hogy a feladatot magam keszitettem, es ha barmilyen segitseget igenybe vettem vagy
// mas szellemi termeket felhasznaltam, akkor a forrast es az atvett reszt kommentekben egyertelmuen jeloltem.
// A forrasmegjeloles kotelme vonatkozik az eloadas foliakat es a targy oktatoi, illetve a
// grafhazi doktor tanacsait kiveve barmilyen csatornan (szoban, irasban, Interneten, stb.) erkezo minden egyeb
// informaciora (keplet, program, algoritmus, stb.). Kijelentem, hogy a forrasmegjelolessel atvett reszeket is ertem,
// azok helyessegere matematikai bizonyitast tudok adni. Tisztaban vagyok azzal, hogy az atvett reszek nem szamitanak
// a sajat kontribucioba, igy a feladat elfogadasarol a tobbi resz mennyisege es minosege alapjan szuletik dontes.
// Tudomasul veszem, hogy a forrasmegjeloles kotelmenek megsertese eseten a hazifeladatra adhato pontokat
// negativ elojellel szamoljak el es ezzel parhuzamosan eljaras is indul velem szemben.
//=============================================================================================

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>

#if defined(__APPLE__)
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__)
#include <windows.h>
#endif
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Innentol modosithatod...

//--------------------------------------------------------
// 3D Vektor
//--------------------------------------------------------
struct Vector {
    float x, y, z;

    Vector() {
        x = y = z = 0;
    }

    Vector(float x0, float y0, float z0 = 0) {
        x = x0;
        y = y0;
        z = z0;
    }

    Vector operator*(float a) {
        return Vector(x * a, y * a, z * a);
    }

    Vector operator/(float a) {
        return Vector(x / a, y / a, z / a);
    }

    Vector operator+(const Vector &v) {
        return Vector(x + v.x, y + v.y, z + v.z);
    }

    Vector operator-(const Vector &v) {
        return Vector(x - v.x, y - v.y, z - v.z);
    }

    float operator*(const Vector &v) {    // dot product
        return (x * v.x + y * v.y + z * v.z);
    }

    Vector operator%(const Vector &v) {    // cross product
        return Vector(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
    }

    float Length() {
        return (float) sqrt(x * x + y * y + z * z);
    }
};

//--------------------------------------------------------
// Spektrum illetve szin
//--------------------------------------------------------
struct Color {
    float r, g, b;

    Color() {
        r = g = b = 0;
    }

    Color(float r0, float g0, float b0) {
        r = r0;
        g = g0;
        b = b0;
    }

    Color operator*(float a) {
        return Color(r * a, g * a, b * a);
    }

    Color operator*(const Color &c) {
        return Color(r * c.r, g * c.g, b * c.b);
    }

    Color operator+(const Color &c) {
        return Color(r + c.r, g + c.g, b + c.b);
    }
};

const int screenWidth = 600;    // alkalmazás ablak felbontása
const int screenHeight = 600;

const size_t maxControlPoints = 10;
const float circleRadius = 2.0;

// színek
const Color BLACK = Color(0.0f, 0.0f, 0.0f);
const Color GREY = Color(0.85f, 0.85f, 0.85f);
const Color TURQUOISE = Color(0.65f, 0.96f, 0.94f);
const Color RED = Color(0.5f, 0.0f, 0.0f);
const Color GREEN = Color(0.0f, 0.5f, 0.0f);
const Color BLUE = Color(0.0f, 0.0f, 0.5f);

long currentTime;
long circulationStartTime;

enum possibleStates {
    ADDING_POINTS, CIRCULATE, CAMERA_MOVING
} currentSate;

struct ControlPoint {
    Vector originalP;
    Vector p;
    float t;
};

ControlPoint cp[maxControlPoints];
size_t cpSize;

class Shape {
protected:
    static const size_t shapeResolution = 200;
    Vector shapePoints[shapeResolution + 1];
    size_t shapePointSize;
    Color color;

public:
    Shape() {
        shapePointSize = 0;
    }

    virtual void computeShape() {
    }

    virtual void drawShape() {
        glColor3f(color.r, color.g, color.b);
        glBegin(GL_LINE_STRIP);
        for (size_t i = 0; i < shapePointSize; i++)
            glVertex2f(shapePoints[i].x, shapePoints[i].y);
        glEnd();
    }
};

class ConvexHull : public Shape {
    bool isClockWise(Vector *p1, Vector *p2, Vector *p3) {
        Vector side1 = *p2 - *p1;
        Vector side2 = *p3 - *p1;
        Vector meroleges = side1 % side2;
        float direction = meroleges.z;
        return (direction <= 0.0);
    }

    static int compareVectorsByX(const void *v1, const void *v2) {
        float x1 = (*(Vector **) v1)->x;
        float x2 = (*(Vector **) v2)->x;
        if (x1 != x2)
            return x1 > x2;

        else {
            float y1 = (*(Vector **) v1)->y;
            float y2 = (*(Vector **) v2)->y;
            return y1 > y2;
        }
    }

    // Az algoritmust A. M. Andrew publikálta 1979-ben. Az implementációt én csináltam.
    void monotoneChain() {
        //  Sort the points of P by x-coordinate (in case of a tie, sort by y-coordinate).
        Vector *sortedCp[cpSize];
        for (size_t i = 0; i < cpSize; i++)
            sortedCp[i] = &cp[i].p;
        qsort(sortedCp, cpSize, sizeof(Vector *), compareVectorsByX);

        //  Initialize U and L as empty lists.
        //  The lists will hold the vertices of upper and lower hulls respectively.
        Vector *lowerHull[cpSize];
        size_t lowerHullSize = 0;
        Vector *upperHull[cpSize];
        size_t upperHullSize = 0;


        //  for i = 1, 2, ..., n:
        //  while L contains at least two points and the sequence of last two points
        //  of L and the point P[i] does not make a counter-clockwise turn:
        //  remove the last point from L
        //  append P[i] to L
        for (size_t i = 0; i < cpSize; i++) {
            while (lowerHullSize >= 2
                    && isClockWise(lowerHull[lowerHullSize - 2], lowerHull[lowerHullSize - 1], sortedCp[i]))
                lowerHullSize--;
            lowerHull[lowerHullSize++] = sortedCp[i];
        }


        //  for i = n, n-1, ..., 1:
        //  while U contains at least two points and the sequence of last two points
        //  of U and the point P[i] does not make a counter-clockwise turn:
        //  remove the last point from U
        //  append P[i] to U
        for (int i = (int) (cpSize - 1); i >= 0; i--) {
            while (upperHullSize >= 2
                    && isClockWise(upperHull[upperHullSize - 2], upperHull[upperHullSize - 1], sortedCp[i]))
                upperHullSize--;
            upperHull[upperHullSize++] = sortedCp[i];
        }

        //  Remove the last point of each list (it's the same as the first point of the other list).
        lowerHullSize--;
        upperHullSize--;

        //  Concatenate L and U to obtain the convex hull of P.
        //  Points in the result will be listed in counter-clockwise order.
        shapePointSize = 0;
        for (size_t i = 0; i < lowerHullSize; i++)
            shapePoints[shapePointSize++] = *lowerHull[i];

        for (size_t i = 0; i < upperHullSize; i++)
            shapePoints[shapePointSize++] = *upperHull[i];
    }

public:
    ConvexHull() : Shape() {
        color = TURQUOISE;
    }

    void computeShape() {
        monotoneChain();
    }

    void drawShape() {
        glColor3f(color.r, color.g, color.b);
        glBegin(GL_TRIANGLE_FAN);
        for (size_t i = 0; i < shapePointSize; i++)
            glVertex2f(shapePoints[i].x, shapePoints[i].y);
        glEnd();
    }

};

class CatmullRomSpline : public Shape {
    unsigned pointsBetweenControlPoints;
    Vector v[maxControlPoints];
    Vector startV;
    Vector endV;

    Vector a0(size_t prev) {
        return cp[prev].p;
    }

    Vector a1(size_t prev) {
        return v[prev];
    }

    Vector a2(size_t prev) {
        size_t i = prev;
        Vector p0 = cp[i].p;
        Vector p1 = cp[i + 1].p;
        float t0 = cp[i].t;
        float t1 = cp[i + 1].t;
        Vector tag1 = (p1 - p0) * 3
                / pow(t1 - t0, 2);
        Vector tag2 = (v[i + 1] + v[i] * 2)
                / (t1 - t0);

        return tag1 - tag2;
    }

    Vector a3(size_t prev) {
        size_t i = prev;
        Vector p0 = cp[i].p;
        Vector p1 = cp[i + 1].p;
        float t0 = cp[i].t;
        float t1 = cp[i + 1].t;
        Vector tag1 = (p0 - p1) * 2
                / pow(t1 - t0, 3);
        Vector tag2 = (v[i + 1] + v[i])
                / pow(t1 - t0, 2);

        return tag1 + tag2;
    }

public:
    CatmullRomSpline() : Shape() {
        startV = Vector(0.00001, 0.00001, 0.0);
        endV = Vector(0.00001, 0.00001, 0.0);
        pointsBetweenControlPoints = shapeResolution / maxControlPoints;
        color = GREEN;
    }

    void computeV() {
        v[0] = startV;
        v[cpSize - 1] = endV;
        for (size_t i = 1; i < cpSize - 1; i++) {
            Vector p0 = cp[i].p;
            Vector pp1 = cp[i + 1].p;
            Vector pm1 = cp[i - 1].p;
            float t0 = cp[i].t;
            float tp1 = cp[i + 1].t;
            float tm1 = cp[i - 1].t;

            Vector tag1 = (pp1 - p0) / (tp1 - t0);
            Vector tag2 = (p0 - pm1) / (t0 - tm1);
            v[i] = tag1 + tag2;
        }
    }

    Vector getPos(float t, size_t prevIndex) {
        size_t i = prevIndex;
        Vector ai0 = a0(i);
        Vector ai1 = a1(i);
        Vector ai2 = a2(i);
        Vector ai3 = a3(i);

        float t0 = cp[i].t;

        return ai3 * pow(t - t0, 3)
                + ai2 * pow(t - t0, 2)
                + ai1 * (t - t0)
                + ai0;
    }

    void computeShape() {
        unsigned points = pointsBetweenControlPoints;
        shapePointSize = 0;

        computeV();
        for (size_t i = 0; i < cpSize - 1; i++)
            for (size_t j = i * points; j < (i + 1) * points; j++) {
                float t = cp[i].t + (
                        ((cp[i + 1].t - cp[i].t) / (float) points) * (j - (i * (float) points))
                );
                shapePoints[j] = getPos(t, i);
                shapePointSize = j + 1;
            }
    }
};

class BezierCurve : public Shape {
    float B (int i, float t) {
        int n = (int) (cpSize - 1);
        float choose = 1;
        for (int j = 1; j <= i; j++)
            choose *=(float) (n - j + 1) / j;
        return (float) (choose * pow(t, i) * pow(1 - t, n - i));
    }

    Vector r (float t) {
        Vector rr(0, 0);
        for (int i = 0; i < cpSize; i++)
            rr = rr + cp[i].p * B(i, t);
        return rr;
    }

public:
    BezierCurve() : Shape() {
        color = RED;
    };

    void computeShape() {
        float range = cp[cpSize - 1].t - cp[0].t;
        shapePointSize = 0;
        for (float t = cp[0].t; t < cp[cpSize - 1].t; t += range / (float) shapeResolution) {
            float t1 = (t - cp[0].t) / range;
            shapePoints[shapePointSize++] = r(t1);
        }
    }
};

class CatmullClarkCurve : public Shape {
    Vector halfSegment(Vector p1, Vector p2) {
        return (p1 + p2) / 2.0;
    }

    Vector centroidOfTriangle(Vector p1, Vector p2, Vector p3, float weightP1 = 0.3333, float weightP2 = 0.3333, float weightP3 = 0.3333) {
        return p1 * weightP1 + p2 * weightP2 + p3 * weightP3;
    }

    void copyArray(Vector from[], size_t fromSize, Vector to[], size_t *toSize) {
        for (size_t i = 0; i < fromSize; i++)
            to[i] = from[i];
        *toSize = fromSize;
    }

    void computeSegmentHalves(Vector *oldCurve, size_t oldCurveSize, Vector *halves) {
        for (size_t i = 1; i < oldCurveSize; i++)
            halves[i - 1] = halfSegment(oldCurve[i - 1], oldCurve[i]);
    }

    void computeCentroids(Vector oldCurve[], size_t oldCurveSize, Vector halves[], Vector centroids[]) {
        size_t halvesSize = oldCurveSize - 1;
        // az 1. és utolsó pont nem változik
        centroids[0] = oldCurve[0];
        centroids[oldCurveSize - 1] = oldCurve[oldCurveSize - 1];
        //  az 1. háromszög: halves[0], oldCurve[1], halves[1]
        // utolsó háromszög: halves[oldCurveSize - 3], oldCurve[oldCurveSize - 2], halves[oldCurveSize - 2]
        for (size_t i = 1; i < halvesSize; i++)
            centroids[i] = centroidOfTriangle(halves[i - 1], oldCurve[i], halves[i], 0.25, 0.5, 0.25);
    }

    // két tömböt össze-merge-el, egyszer az elsőből vesz elemet, aztán a másikból
    void mergeAlternating(Vector arr1[], size_t arr1Size, Vector arr2[], size_t arr2Size, Vector newArr[]) {
        // arr1: centroids, oldCurveSize
        // arr2: halves, halvesSize
        // 0 -> newCurveSize - 1 = oldCurveSize  * 2 - 2
        for (size_t i = 0; i < arr1Size; i++)
            newArr[2 * i] = arr1[i];
        // 1 -> newCurveSize - 2 = oldCurveSize * 2 - 3
        for (size_t i = 0; i < arr2Size; i++)
            newArr[2 * i + 1] = arr2[i];
    }

public:
    CatmullClarkCurve() : Shape() {
        color = BLUE;
    };

    void computeShape() {
        Vector oldCurve[shapeResolution + 1];
        Vector newCurve[shapeResolution + 1];
        size_t newCurveSize;
        size_t oldCurveSize;

        // a new curve-öt a kontrollpontokra inicializálja
        newCurveSize = cpSize;
        for (size_t i = 0; i < cpSize; i++)
            newCurve[i] = cp[i].p;

        while (newCurveSize < shapeResolution / 2) {
            // átmásolja newCurve-öt oldCurve-be
            copyArray(newCurve, newCurveSize, oldCurve, &oldCurveSize);

            // kiszámítja az old Curve szakaszainak felezőpontjait, 1-el kevesebb felezőpont van, mint pont.
            size_t halvesSize = oldCurveSize - 1;
            Vector halves[halvesSize];
            computeSegmentHalves(oldCurve, oldCurveSize, halves);

            // kiszámítja, hová kell eltolni az old Curve pontjait, hogy a háromszög súlypontjába kerüljenek
            Vector centroids[oldCurveSize];
            computeCentroids(oldCurve, oldCurveSize, halves, centroids);

            // merge: centroids + halves
            mergeAlternating(centroids, oldCurveSize, halves, halvesSize, newCurve);
            newCurveSize = oldCurveSize + halvesSize;
        }

        // lecseréli shapePoints-ot newCurve-re
        copyArray(newCurve, newCurveSize, shapePoints, &shapePointSize);
    }
};

class Camera {
    static const float windowWidth = 58;
    static const float windowHeight = 68;
    Vector bottomLeft, topRight, originalBottomLeft, originalTopRight;
public:
    Camera() {
        originalBottomLeft = bottomLeft = Vector(21.0f, 16.0f);
        originalTopRight = topRight = bottomLeft + Vector(windowWidth, windowHeight);
    }

    float left() {
        return bottomLeft.x;
    }

    float right() {
        return topRight.x;
    }

    float top() {
        return topRight.y;
    }

    float bottom() {
        return bottomLeft.y;
    }

    void move(Vector change) {
        bottomLeft = originalBottomLeft + change;
        topRight = originalTopRight + change;
    }

    Vector getWorldPos(int x, int y, int screenWidth, int screenHeight) {
        float x01 = (float) x / (float) screenWidth;
        float yInverse = screenHeight - y;
        float y01 = yInverse / (float) screenHeight;
        float xResized = x01 * windowWidth;
        float yResized = y01 * windowHeight;
        Vector resizedPos(xResized, yResized, 0);
        return resizedPos + bottomLeft;
    }

} camera;

void addControlPoint(Vector p, float t) {
    if (cpSize < maxControlPoints) {
        cp[cpSize].p = p;
        cp[cpSize].originalP = p;
        cp[cpSize].t = t;
        cpSize++;
    }
}

ControlPoint *getControlPointAtPos(Vector p) {
    for (size_t i = 0; i < cpSize; i++) {
        Vector distance = p - cp[i].p;
        if (distance.Length() < circleRadius)
            return &cp[i];
    }
    return NULL;
}

ControlPoint *followedControlPoint;
Vector followStartPoint;
ConvexHull convexHull;
CatmullRomSpline crSpline;
BezierCurve bezier;
CatmullClarkCurve clark;
Shape *shapes[4] = {&convexHull, &crSpline, &bezier, &clark};

// Inicializacio, a program futasanak kezdeten, az OpenGL kontextus letrehozasa utan hivodik meg (ld. main() fv.)
void onInitialization() {
    glViewport(0, 0, screenWidth, screenHeight);
    currentSate = ADDING_POINTS;
    cpSize = 0;
}

void drawCircle(Vector Center, float radius) {
    int triangles = 30;
    float Pi = 3.14159;
    float delta = 2 * Pi / triangles;
    glBegin(GL_TRIANGLE_FAN);
    glVertex2f(Center.x, Center.y);
    for (int i = 0; i <= triangles; i++)
        glVertex2f(
                Center.x + (radius * (float) cos(i * delta)),
                Center.y + (radius * (float) sin(i * delta))
        );
    glEnd();
}

void moveControlPoints(float ts, long circulationStartTime) {
    double period = 5000.0;
    double radius = 5;
    double Pi = 3.14159;
    if (ts < circulationStartTime)
        ts = circulationStartTime;
    // azért adjuk hozzá, hogy a kör tetejéről induljon a periódus, a periódus 1/4-edénél van a kör tetején.
    ts += 0.25 * period;
    double periodsDone = floor((ts - circulationStartTime) / period);
    double timePastInPeriod = ts - circulationStartTime - (periodsDone * period);
    double periodPart = timePastInPeriod / period;
    double angle = 2 * Pi * periodPart;

    for (size_t i = 0; i < cpSize; i += 2) {
        float newX = (float) (cp[i].originalP.x - radius * cos(angle));
        float newY = (float) (cp[i].originalP.y + radius * sin(angle) - radius);
        cp[i].p = Vector(newX, newY);
    }

    for (size_t i = 1; i < cpSize; i += 2) {
        float newX = (float) (cp[i].originalP.x + radius * cos(angle));
        float newY = (float) (cp[i].originalP.y + radius * sin(angle) - radius);
        cp[i].p = Vector(newX, newY);
    }
}

void SimulateWorld(float tstart, float tend) {
    float dt = 10.0f;
    if (currentSate == CIRCULATE || currentSate == CAMERA_MOVING)
        for (float ts = tstart; ts < tend; ts += dt) {
            moveControlPoints(ts, circulationStartTime);
            if (cpSize >= 2)
                for (size_t i = 0; i < 4; i++)
                    shapes[i]->computeShape();
        }

    if (currentSate == CAMERA_MOVING)
        camera.move(followedControlPoint->p - followStartPoint);
}

// Rajzolas, ha az alkalmazas ablak ervenytelenne valik, akkor ez a fuggveny hivodik meg
void onDisplay() {
    glClearColor(GREY.r, GREY.g, GREY.b, 1.0f);     // torlesi szin beallitasa
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // kepernyo torles

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(camera.left(), camera.right(), camera.bottom(), camera.top());


    if (cpSize >= 2) {
        for (size_t i = 0; i < 4; i++)
            shapes[i]->drawShape();

        glColor3f(BLACK.r, BLACK.g, BLACK.b);
        for (size_t i = 0; i < cpSize; i++)
            drawCircle(cp[i].p, circleRadius);
    }

    else {
        glColor3f(BLACK.r, BLACK.g, BLACK.b);
        for (size_t i = 0; i < cpSize; i++)
            drawCircle(cp[i].p, circleRadius);
    }

    glutSwapBuffers();                    // Buffercsere: rajzolas vege

}

// Billentyuzet esemenyeket lekezelo fuggveny (lenyomas)
void onKeyboard(unsigned char key, int x, int y) {
    // space-re elindítjuk a keringést, de csak ha még a keringés előtti állapotban vagyunk
    // amúgy nem törődünk a Space-el
    if (key == ' ' && currentSate == ADDING_POINTS) {
        circulationStartTime = glutGet(GLUT_ELAPSED_TIME);
        currentSate = CIRCULATE;
    }
}

// Billentyuzet esemenyeket lekezelo fuggveny (felengedes)
void onKeyboardUp(unsigned char key, int x, int y) {

}

// Eger esemenyeket lekezelo fuggveny
void onMouse(int button, int state, int x, int y) {
    if (button == GLUT_LEFT_BUTTON && state == GLUT_UP && currentSate == ADDING_POINTS) {
        long clickTime = glutGet(GLUT_ELAPSED_TIME);

        // felvesz egy új kontrollpontot és újraszámolja a görbét, de csak a pontfelvétel állapotában, ha még nem értük el a maxot
        Vector pos = camera.getWorldPos(x, y, screenWidth, screenHeight);
        if (cpSize < maxControlPoints) {
            addControlPoint(pos, clickTime / 1000.0f);
            if (cpSize >= 2)
                for (size_t i = 0; i < 4; i++)
                    shapes[i]->computeShape();
            glutPostRedisplay();
        }
    }

    if (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN) {
        if (currentSate == CIRCULATE || currentSate == CAMERA_MOVING) {
            Vector pos = camera.getWorldPos(x, y, screenWidth, screenHeight);
            ControlPoint *clickedControlPoint = getControlPointAtPos(pos);
            if (clickedControlPoint != NULL) {
                followedControlPoint = clickedControlPoint;
                followStartPoint = clickedControlPoint->p;
                camera.move(followedControlPoint->p - followStartPoint);
                currentSate = CAMERA_MOVING;
            }
        }
        glutPostRedisplay();
    }
}

// Eger mozgast lekezelo fuggveny
void onMouseMotion(int x, int y) {

}


// `Idle' esemenykezelo, jelzi, hogy az ido telik, az Idle esemenyek frekvenciajara csak a 0 a garantalt minimalis ertek
void onIdle() {
    float old_time = currentTime;
    currentTime = glutGet(GLUT_ELAPSED_TIME);        // program inditasa ota eltelt ido

    SimulateWorld(old_time, currentTime);
    glutPostRedisplay();
}


// ...Idaig modosithatod
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// A C++ program belepesi pontja, a main fuggvenyt mar nem szabad bantani
int main(int argc, char **argv) {
    glutInit(&argc, argv); 				// GLUT inicializalasa
    glutInitWindowSize(600, 600);			// Alkalmazas ablak kezdeti merete 600x600 pixel
    glutInitWindowPosition(100, 100);			// Az elozo alkalmazas ablakhoz kepest hol tunik fel
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);	// 8 bites R,G,B,A + dupla buffer + melyseg buffer

    glutCreateWindow("Grafika hazi feladat");		// Alkalmazas ablak megszuletik es megjelenik a kepernyon

    glMatrixMode(GL_MODELVIEW);				// A MODELVIEW transzformaciot egysegmatrixra inicializaljuk
    glLoadIdentity();
    glMatrixMode(GL_PROJECTION);			// A PROJECTION transzformaciot egysegmatrixra inicializaljuk
    glLoadIdentity();

    onInitialization();					// Az altalad irt inicializalast lefuttatjuk

    glutDisplayFunc(onDisplay);				// Esemenykezelok regisztralasa
    glutMouseFunc(onMouse);
    glutIdleFunc(onIdle);
    glutKeyboardFunc(onKeyboard);
    glutKeyboardUpFunc(onKeyboardUp);
    glutMotionFunc(onMouseMotion);

    glutMainLoop();					// Esemenykezelo hurok

    return 0;
}