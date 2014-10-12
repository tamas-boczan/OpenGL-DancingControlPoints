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
#include <stdio.h>

#endif


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Innentol modosithatod...

//--------------------------------------------------------
// 3D Vektor
//--------------------------------------------------------
struct Vector {
    float x, y, z;

    Vector( ) {
        x = y = z = 0;
    }
    Vector(float x0, float y0, float z0 = 0) {
        x = x0; y = y0; z = z0;
    }
    Vector operator*(float a) {
        return Vector(x * a, y * a, z * a);
    }
    Vector operator/(float a) {
        return Vector(x / a, y / a, z / a);
    }
    Vector operator+(const Vector& v) {
        return Vector(x + v.x, y + v.y, z + v.z);
    }
    Vector operator-(const Vector& v) {
        return Vector(x - v.x, y - v.y, z - v.z);
    }
    float operator*(const Vector& v) { 	// dot product
        return (x * v.x + y * v.y + z * v.z);
    }
    Vector operator%(const Vector& v) { 	// cross product
        return Vector(y*v.z-z*v.y, z*v.x - x*v.z, x*v.y - y*v.x);
    }
    float Length() {
        return sqrt(x * x + y * y + z * z);
    }

    Vector normalized() {
        return Vector(x / Length(), y / Length(), z / Length());
    }
};

//--------------------------------------------------------
// Spektrum illetve szin
//--------------------------------------------------------
struct Color {
    float r, g, b;

    Color( ) {
        r = g = b = 0;
    }
    Color(float r0, float g0, float b0) {
        r = r0; g = g0; b = b0;
    }
    Color operator*(float a) {
        return Color(r * a, g * a, b * a);
    }
    Color operator*(const Color& c) {
        return Color(r * c.r, g * c.g, b * c.b);
    }
    Color operator+(const Color& c) {
        return Color(r + c.r, g + c.g, b + c.b);
    }
};

const int screenWidth = 600;	// alkalmazás ablak felbontása
const int screenHeight = 600;

const size_t maxControlPoints = 10;
const size_t curveResolution = 200;

// színek
const Color BLACK = Color(0.0f, 0.0f, 0.0f);
const Color GREY = Color(0.85f, 0.85f, 0.85f);
// Türkíz szín forrása: www.opengl.org/discussion_boards/showthread.php/132502-Color-tables
const Color TURQUOISE = Color(0.678f, 0.918f, 0.918f);
const Color RED = Color(1.0f, 0.0f, 0.0f);
const Color GREEN = Color(0.0f, 0.6f, 0.1f);
const Color BLUE = Color(0.0f, 0.0f, 1.0f);


long currentTime;
long clickTime;
long circulationStartTime;
unsigned int followedControlPoint;

enum possibleStates{ADDING_POINTS, CIRCULATE, CAMERA_MOVING} currentSate;

struct ControlPoint
{
    Vector p;
    float t;
};

ControlPoint allControlPoints[maxControlPoints];
int allControlPointNumber;

class ControlPointList
{
    ControlPoint* controlPoints[maxControlPoints];
    size_t size;

public:

    ControlPointList(size_t size = 0) : size(size) {

    }

    ControlPointList (const ControlPointList &c2)
    {
        size = c2.size;
        for (int i = 0; i < size; i++) {
            controlPoints[i] = c2.controlPoints[i];
        }
    }

    ControlPoint * getControlPoint(unsigned i) {
        return controlPoints[i];
    }

    void setControlPoint(unsigned i, ControlPoint * cp){
        controlPoints[i] = cp;
    }

    Vector getP(unsigned i)
    {
        return controlPoints[i]-> p;
    }

    void setP(unsigned i, Vector Position)
    {
        controlPoints[i] -> p = Position;
    }

    float getT(unsigned i)
    {
        return controlPoints[i] -> t;
    }

    void setT(unsigned i, float weight)
    {
        controlPoints[i] -> t = weight;
    }

    size_t getSize() {
        return size;
    }

    void add(ControlPoint * cp) {
        if (size < maxControlPoints) {
            controlPoints[size] = cp;
            size++;
        }
    }

    void add(Vector _p, float _t) {
        if (size < maxControlPoints) {
            allControlPoints[allControlPointNumber].p = _p;
            allControlPoints[allControlPointNumber].t = _t;
            add(allControlPoints + allControlPointNumber);
            allControlPointNumber++;
        }

    }

    void removeLastReference() {
        size--;
    }

    void clear() {
        size = 0;
    }
};

ControlPointList currentControlPoints;
ControlPointList originalControlPoints;
ControlPointList convexHull;

/*
int findLowestControlPoint (ControlPointList * currentControlPoints){
    int lowestIndex;
    float lowestY = 9999999999999;
    for (unsigned i = 0; i < currentControlPoints->getSize(); i++){
        float currentY = currentControlPoints->getControlPoint(i).p.y;
        if (currentY < lowestY){
            lowestY = currentY;
            lowestIndex = i;
        }
    }
    return lowestIndex;
}
*/

bool isCounterClockWise (Vector p1, Vector p2, Vector p3) {
    Vector side1 = p2 - p1;
    Vector side2 = p3 - p1;
    Vector meroleges = side1 % side2;
    float direction = meroleges.z;
    return (direction <= 0.0);
}

bool isLeftFrom (Vector v1, Vector v2){
    if (v1.x < v2.x)
        return true;
    if (v1.x > v2.x)
        return false;
    return (v1.y < v2.y);
}

void sortControlPointsByX(ControlPointList * cpList){
    size_t n = cpList->getSize();
    for (size_t c = 0 ; c < ( n - 1 ); c++) {
        for (size_t d = 0; d < n - c - 1; d++) {
            if (isLeftFrom(cpList->getP(d + 1), cpList->getP(d))) {
                ControlPoint *temp = cpList->getControlPoint(d);
                cpList->setControlPoint(d, cpList->getControlPoint(d + 1));
                cpList->setControlPoint(d + 1, temp);
            }
        }
    }
}


// Returns a list of points on the convex hull in counter-clockwise order.
// Note: the last point in the returned list is the same as the first one.
ControlPointList findConvexHull (ControlPointList currentControlPoints){
    ControlPointList orderedCp(currentControlPoints);
    sortControlPointsByX(&orderedCp);

    ControlPointList lowerHull;
    for (int i = 0; i < orderedCp.getSize(); i++){
        size_t hullSize = lowerHull.getSize();
        while (hullSize >= 2 && !isCounterClockWise(lowerHull.getP(hullSize - 2), lowerHull.getP(hullSize - 1), orderedCp.getP(i))) {
            lowerHull.removeLastReference();
            hullSize = lowerHull.getSize();
        }
        lowerHull.add(orderedCp.getControlPoint(i));
    }

    /*
    for i = n, n-1, ..., 1:
    while U contains at least two points and the sequence of last two points
            of U and the point P[i] does not make a counter-clockwise turn:
        remove the last point from U
    append P[i] to U
     */
    ControlPointList upperHull;
    for (int i = orderedCp.getSize() - 1; i > 0; i--){
        size_t hullSize = upperHull.getSize();
        while (hullSize >= 2 && !isCounterClockWise(upperHull.getP(hullSize - 2), upperHull.getP(hullSize - 1), orderedCp.getP(i))) {
            upperHull.removeLastReference();
            hullSize = upperHull.getSize();
        }
        upperHull.add(orderedCp.getControlPoint(i));
    }

    /*
    Remove the last point of each list (it's the same as the first point of the other list).
    Concatenate L and U to obtain the convex hull of P.
    Points in the result will be listed in counter-clockwise order.
     */
    if (lowerHull.getControlPoint(lowerHull.getSize() - 1) ==
            upperHull.getControlPoint(0))
        lowerHull.removeLastReference();

    if (upperHull.getControlPoint(upperHull.getSize() - 1) ==
            lowerHull.getControlPoint(0))
        upperHull.removeLastReference();

    convexHull.clear();
    for (int i = 0; i < lowerHull.getSize(); i++)
        convexHull.add(lowerHull.getControlPoint(i));
    for (int i = 0; i < upperHull.getSize(); i++)
        convexHull.add(upperHull.getControlPoint(i));
}


/*
int* findConvexHull (ControlPointList * currentControlPoints){
    bool usedCP[currentControlPoints->getSize()];
    for (int i = 0; i < currentControlPoints->getSize(); i++)
        usedCP = false;

    // legalsó kontrollpont
    int nextIndex = findLowestControlPoint(currentControlPoints);
    convexHullIndexes[0] = nextIndex;
    convexHullSize = 1;
    usedCP[nextIndex] = true;

    //egyik oldalról záró kontrollpontok keresése
    bool found = false;
    int lastIndex = convexHullIndexes[convexHullSize - 1];
    for (int i = 0; i < currentControlPoints->getSize(); i++)
        if (!usedCP[i]){
            // last - i egyenes
            Vector lineV = currentControlPoints->getP(lastIndex) - currentControlPoints->getP(i);
            Vector lineP = currentControlPoints->getP(lastIndex);
            // ellenőrzi, hogy minden pont vagy balra vagy jobbra van-e ettől az egyenestől
        }
}
*/

class Curve
{
protected:
    ControlPointList *cp;
    Vector curvePoints[curveResolution + 1];
    size_t curvePointSize;

public:
    Curve (ControlPointList * controlPoints){
        cp = controlPoints;
        curvePointSize = 0;
    }

    Vector getCurvePoint(unsigned i){
        return curvePoints[i];
    }

    size_t getCurvePointSize(){
        return curvePointSize;
    }

    void drawCurve(Color color) {
        glColor3f(color.r, color.g, color.b);
        glBegin(GL_LINE_STRIP);
        for (unsigned i = 0; i < curvePointSize; i++)
            glVertex2f(curvePoints[i].x, curvePoints[i].y);
        glEnd();
    }
};

class CRSpline: public Curve
{
    unsigned pointsBetweenControlPoints;
    //sebesség
    Vector sebesseg[maxControlPoints];
    Vector kezdosebesseg;
    Vector vegsebesseg;

    Vector GetAi2(int elozo)
    {
        int i = elozo;
        Vector p0 = cp->getP(i);
        Vector p1 = cp->getP(i+1);
        float t0 = cp->getT(i);
        float t1 = cp->getT(i + 1);
        Vector tag1 = (p1 - p0) * 3
                / pow(t1 - t0, 2);
        Vector tag2 = (sebesseg[i + 1] + sebesseg[i] * 2)
                / (t1 - t0);

        return tag1 - tag2;
    }

    Vector GetAi3(int elozo)
    {
        int i = elozo;
        Vector p0 = cp->getP(i);
        Vector p1 = cp->getP(i+1);
        float t0 = cp->getT(i);
        float t1 = cp->getT(i + 1);
        Vector tag1 = (p0 - p1) * 2
                / pow(t1 - t0, 3);
        Vector tag2 = (sebesseg[i + 1] + sebesseg[i])
                / pow(t1 - t0, 2);

        return tag1 + tag2;
    }

public:
    CRSpline(ControlPointList * controlPointList)
            : Curve(controlPointList)    // Call the superclass constructor in the subclass' initialization list.
    {
        kezdosebesseg = Vector(0.00001, 0.00001, 0.0);
        vegsebesseg = Vector(0.00001, 0.00001, 0.0);
        pointsBetweenControlPoints = curveResolution / maxControlPoints;
    }

    void ComputeV()
    {
        sebesseg[0] = kezdosebesseg;
        sebesseg[cp->getSize() - 1] = vegsebesseg;
        for (int i = 1; i < cp->getSize() - 1; i++)
        {
            Vector p0 = cp->getP(i);
            Vector pp1 = cp->getP(i + 1);
            Vector pm1 = cp->getP(i - 1);
            float t0 = cp->getT(i);
            float tp1 = cp->getT(i + 1);
            float tm1 = cp->getT(i - 1);

            Vector tag1 = (pp1 - p0) / (tp1 - t0);
            Vector tag2 = (p0 - pm1) / (t0 - tm1);
            sebesseg[i] = tag1 + tag2;
        }
    }

    Vector GetPos(float t, int elozo_kontrollpont_szama)
    {
        int i = elozo_kontrollpont_szama;
        Vector ai0 = cp->getP(i);
        Vector ai1 = sebesseg[i];
        Vector ai2 = GetAi2(i);
        Vector ai3 = GetAi3(i);

        float t0 = cp->getT(i);

        return ai3 * pow(t - t0, 3)
                + ai2 * pow(t - t0, 2)
                + ai1 * (t - t0)
                + ai0;
    }

    Vector GetPos(float t)
    {
        int index=0;
        for (int i = 0; cp->getT(i) < t; i++)
            index = i;
        return GetPos(t, index);
    }

    Vector derivaltPos(float t, int elozo_kontrollpont_szama)
    {
        int i = elozo_kontrollpont_szama;
        Vector ai1 = sebesseg[i];
        Vector ai2 = GetAi2(i);
        Vector ai3 = GetAi3(i);

        return ai3 * pow(t - cp->getT(i), 2) * 3
                + ai2 * (t - cp->getT(i)) * 2
                + ai1;
    }

    Vector Tangencialis(float t, int index)
    {
        Vector T = derivaltPos(t, index);
        return T.normalized();
    }

    Vector Binormalis(float t, int index)
    {
        Vector B = Tangencialis(t, index) % Vector(1, 0, 0);
        return B.normalized();
    }

    Vector Normal(float t, int index)
    {
        Vector N = Tangencialis(t, index) % Binormalis(t, index);
        return N.normalized();
    }

    void computeCurve() {
        size_t cpSize = cp->getSize();
        unsigned points = pointsBetweenControlPoints;
        curvePointSize = 0;

        ComputeV();
        for (unsigned i = 0; i < cpSize -1; i++)
            for (unsigned j = i * points; j < (i+1) * points; j++)
            {
                float t = cp->getT(i) + (
                        ((cp->getT(i+1) - cp->getT(i)) / (float)points) * (j - (i * (float)points))
                );
                curvePoints[j] = GetPos(t, i);
                curvePointSize = j + 1;
            }
    }
} crSpline(&currentControlPoints);


class BezierCurve : public Curve
{
    // t = [0-1] tartomany kozott, i=m=kpMax-1
    Vector CountBezierPos(float t, int i, int m)
    {
        if (m == 0)
            return cp->getP(i);

        return CountBezierPos(t, i - 1, m - 1) * t +
                CountBezierPos(t, i, m - 1) * (1.0 - t);
    }

public:
    BezierCurve (ControlPointList * controlPointList)
            : Curve(controlPointList){
        curvePointSize = 0;
    };

    Vector GetPos(float t)
    {
        return CountBezierPos(t, cp->getSize() - 1, cp->getSize() - 1);
    }

    void computeCurve() {
        float firstT = cp->getT(0);
        float lastT = cp->getT(cp->getSize() - 1);
        float range = lastT - firstT;
        int i = 0;
        for (float t = firstT; t < lastT; t += range / (float) curveResolution) {
            float t1 = (t - firstT) / range;
            curvePoints[i] = GetPos(t1);
            curvePointSize = i;
            i++;
        }
    }

} bezier(&currentControlPoints);

class Camera
{
    Vector bottomLeft, topRight;
    Vector originalBottomLeft, originalTopRight;
public:
    Camera()
    {

        originalBottomLeft = bottomLeft = Vector(21.0f, 16.0f);
        originalTopRight = topRight = bottomLeft + Vector(58.0f, 68.0f);
    }

    float width() {
        return topRight.x - bottomLeft.x;
    }

    float height() {
        return topRight.y - bottomLeft.y;
    }

    float left(){
        return bottomLeft.x;
    }

    float right(){
        return topRight.x;
    }

    float top(){
        return topRight.y;
    }

    float bottom(){
        return bottomLeft.y;
    }

    void moveXTo(float movement)
    {
        bottomLeft.x = originalBottomLeft.x + movement;
        topRight.x = originalTopRight.x + movement;
    }

    void moveYTo(float movement)
    {
        bottomLeft.y = originalBottomLeft.y + movement;
        topRight.y = originalTopRight.y + movement;
    }

    Vector getWorldPos (int x, int y, int screenWidth, int screenHeight)
    {
        float x01 = (float) x / (float) screenWidth;
        float yInverse = screenHeight - y;
        float y01 = (float) yInverse / (float) screenHeight;
        float xResized = x01 * Camera::width();
        float yResized = y01 * Camera::height();
        Vector resizedPos(xResized, yResized, 0);
        return resizedPos + bottomLeft;
    }

} camera;

// Inicializacio, a program futasanak kezdeten, az OpenGL kontextus letrehozasa utan hivodik meg (ld. main() fv.)
void onInitialization( )
{
    glViewport(0, 0, screenWidth, screenHeight);
    currentSate = ADDING_POINTS;
    allControlPointNumber = 0;

}

void DrawCircle(Vector Center, float radius)
{
    int triangles = 30;
    float Pi = 3.14159;
    float delta = 2 * Pi / triangles;
    glBegin(GL_TRIANGLE_FAN);
    glVertex2f(Center.x, Center.y);
    for (int i = 0; i <= triangles; i++)
        glVertex2f (
                Center.x + (radius * cos(i *  delta)),
                Center.y + (radius * sin(i * delta))
        );
    glEnd();
}

void moveControlPoints(float ts, long circulationStartTime){
    // float t = (ts/1000.f) - (circleStartTime/1000.0f) + DzsCurve.getT(0);
}

void moveCamera(Vector lookat) {

}

void SimulateWorld(float tstart, float tend)
{
    float dt = 10.0f;
    if (currentSate == CIRCULATE)
        for (float ts = tstart; ts < tend; ts += dt) {
            moveControlPoints(ts, circulationStartTime);

            if (currentControlPoints.getSize() >= 2) {
                crSpline.computeCurve();
                bezier.computeCurve();
                findConvexHull(currentControlPoints);
            }
        }

    if (currentSate == CAMERA_MOVING)
        for (float ts = tstart; ts < tend; ts += dt) {
            moveCamera(currentControlPoints.getP(followedControlPoint));

            if (currentControlPoints.getSize() >= 2) {
                crSpline.computeCurve();
                bezier.computeCurve();
                findConvexHull(currentControlPoints);
            }
        }
}


// Rajzolas, ha az alkalmazas ablak ervenytelenne valik, akkor ez a fuggveny hivodik meg
void onDisplay( )
{
    glClearColor(GREY.r, GREY.g, GREY.b, 1.0f);     // torlesi szin beallitasa
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // kepernyo torles

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(camera.left(), camera.right(), camera.bottom(), camera.top());


    if (currentControlPoints.getSize() >= 2)
    {
        glColor3f(TURQUOISE.r, TURQUOISE.g, TURQUOISE.b);
        glBegin(GL_TRIANGLE_FAN);
        for (unsigned i = 0; i < convexHull.getSize(); i++)
            glVertex2f(convexHull.getP(i).x, convexHull.getP(i).y);
        glEnd();

        crSpline.drawCurve(GREEN);
        bezier.drawCurve(RED);

        glColor3f(BLACK.r, BLACK.g, BLACK.b);
        for (unsigned i = 0; i < currentControlPoints.getSize(); i++)
            DrawCircle(currentControlPoints.getP(i), 2.0f);
    }

    else {
        glColor3f(BLACK.r, BLACK.g, BLACK.b);
        for (unsigned i = 0; i < currentControlPoints.getSize(); i++)
            DrawCircle(currentControlPoints.getP(i), 2.0f);
    }

    glutSwapBuffers();     				// Buffercsere: rajzolas vege

}

// Billentyuzet esemenyeket lekezelo fuggveny (lenyomas)
void onKeyboard(unsigned char key, int x, int y)
{
    // space-re elindítjuk a keringést, de csak ha még a keringés előtti állapotban vagyunk
    // amúgy nem törődünk a Space-el
    if (key == ' ')
    {
        if (currentSate == ADDING_POINTS) {
            circulationStartTime = glutGet(GLUT_ELAPSED_TIME);
            currentSate = CIRCULATE;
        }
    }
}

// Billentyuzet esemenyeket lekezelo fuggveny (felengedes)
void onKeyboardUp(unsigned char key, int x, int y)
{

}

// Eger esemenyeket lekezelo fuggveny
void onMouse(int button, int state, int x, int y)
{
    if (button == GLUT_LEFT_BUTTON && state == GLUT_UP)
    {
        clickTime = glutGet( GLUT_ELAPSED_TIME );

        // felvesz egy új kontrollpontot és újraszámolja a görbét, de csak a pontfelvétel állapotában, ha még nem értük el a maxot
        if (currentSate == ADDING_POINTS){
            Vector pos = camera.getWorldPos(x, y, screenWidth, screenHeight);
            if (currentControlPoints.getSize() < maxControlPoints)
            {
                currentControlPoints.add(pos, clickTime/1000.0f);
                if (currentControlPoints.getSize() >= 2) {
                    crSpline.computeCurve();
                    bezier.computeCurve();
                    findConvexHull(currentControlPoints);
                }

                glutPostRedisplay();
            }
        }
    }

    if (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN)
    {
        // elindítja a kamerát,
        if (currentSate == CIRCULATE || currentSate == CAMERA_MOVING){
            currentSate = CAMERA_MOVING;
            // TODO kitalálni melyik kontrollpontra kattintott
            // followedControlPoint =
        }
        glutPostRedisplay();
    }


}

// Eger mozgast lekezelo fuggveny
void onMouseMotion(int x, int y)
{

}


// `Idle' esemenykezelo, jelzi, hogy az ido telik, az Idle esemenyek frekvenciajara csak a 0 a garantalt minimalis ertek
void onIdle( )
{
    float  old_time = currentTime;
    currentTime = glutGet(GLUT_ELAPSED_TIME);		// program inditasa ota eltelt ido

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
