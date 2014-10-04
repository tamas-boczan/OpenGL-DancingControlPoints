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


class ControlPointList {
    ControlPoint controlPoints [maxControlPoints];
    size_t size;

public:

    ControlPointList(size_t size = 0) : size(size) {
    }

    Vector getP(unsigned i)
    {
        return  controlPoints[i].p;
    }

    void setP(unsigned i, Vector Position)
    {
        controlPoints[i].p = Position;
    }

    float getT(unsigned i)
    {
        return controlPoints[i].t;
    }

    void setT(unsigned i, float weight)
    {
        controlPoints[i].t = weight;
    }

    size_t getSize() {
        return size;
    }

    void add(ControlPoint cp) {
        if (size < maxControlPoints) {
            controlPoints[size] = cp;
            size++;
        }
    }

    void add(Vector p, float t) {
        if (size < maxControlPoints) {
            controlPoints[size].p = p;
            controlPoints[size].t = t;
            size++;
        }
    }
} controlPoints;

class Curve {
protected:
    ControlPointList *cp;
    Vector* curvePoints;
    size_t curvePointSize;

public:
    Curve (ControlPointList * controlPoints){
        cp = controlPoints;
        curvePoints = NULL;
        curvePointSize = 0;
    }

    ~Curve(){
        delete[] curvePoints;
    }

    Vector getCurvePoint(unsigned i){
        return curvePoints[i];
    }

    size_t getCurvePointSize(){
        return curvePointSize;
    }
};

class DzsCurve : public Curve {
    Vector v [maxControlPoints];
    Vector a [maxControlPoints];

public:
    DzsCurve(ControlPointList *controlPoints) : Curve(controlPoints) {
    }

    static const unsigned CurvePtNr = 20;

private:

    Vector a0(int i)
    {
        return cp->getP(i);
    }

    Vector a1(int i)
    {
        return v[i];
    }

    Vector a2(int i)
    {
        return a[i] / 2.0f;
    }

    Vector a3(int i)
    {
        Vector pi = cp -> getP(i);
        Vector pi1 = cp -> getP(i +1);
        float ti = cp->getT(i);
        float ti1 = cp->getT(i + 1);
        Vector vi = v[i];
        Vector vi1 = v[i+1];
        Vector ai = a[i];

        return (
                ((pi1 * 4.0f - pi * 4.0f) / (pow(ti1 - ti, 3)))
              - ((vi1 + vi * 3.0f) / (pow(ti1 - ti, 2)))
                - (ai / (ti1 - ti))
        );
    }

    Vector a4(int i)
    {
        Vector pi = cp -> getP(i);
        Vector pi1 = cp -> getP(i +1);
        float ti = cp->getT(i);
        float ti1 = cp->getT(i + 1);
        Vector vi = v[i];
        Vector vi1 = v[i+1];
        Vector ai = a[i];

        return (
                (((pi1 * (-3.0f)) + (pi * 3.0f)) / (pow((ti1 - ti), 4)))
                        + ((vi1 + (vi * 2.0f)) / (pow((ti1 - ti),3)))
                        + ((ai *0.5f) / (pow((ti1 - ti), 2)))
        );
    }

    void computeV()
    {
        size_t cpSize = cp->getSize();
        v[cpSize -1] = Vector(0.0f, 0.0f, 0.0f);
        for (unsigned i = 1; i < cpSize - 1; i++) {
            Vector pi = cp -> getP(i);
            Vector pim = cp -> getP(i - 1);
            Vector pip = cp -> getP(i + 1);
            float ti = cp->getT(i);
            float tim = cp->getT(i - 1);
            float tip = cp->getT(i + 1);

            v[i] = (
                    ((pi - pim) / (ti - tim)) +
                    ((pip - pi) / (tip - ti))
            ) / 2.0f;
        }
    }

    void computeA() {
        for (unsigned i = 0; i < cp->getSize() -1; i++)
        {
            float ti = cp->getT(i);
            float tip = cp->getT(i + 1);
            a[i + 1] = (a4(i) * 12.0f * (pow((tip - ti),2))) +
                    (a3(i) *  6.0f * (tip - ti)) +
                    (a2(i) *  2.0f);
        }
    }

    size_t findCP(float t)
    {
        for(unsigned i = 1; i < cp->getSize(); i++)
            if (cp->getT(i) > t)
                return i-1;
        return cp->getSize()-1;
    }

public:
    Vector getCurvePos(unsigned previousControlPointIndex, float weight)
    {
        unsigned i = previousControlPointIndex;
        float t = weight;
        float ti = cp->getT(i);

        return	(a4(i) * pow((t-ti),4)) +
                (a3(i) * pow((t-ti),3)) +
                (a2(i) * pow((t-ti),2)) +
                (a1(i) * (t-ti)) +
                a0(i);
    }

    Vector getCurvePos(float weight)
    {
        int previousControlPointIndex = findCP(weight);
        return getCurvePos(previousControlPointIndex, weight);
    }

    Vector getCurvePoint(int i)
    {
        return curvePoints[i];
    }

    void computeCurve()
    {
        computeV();
        computeA();

        size_t cpSize = cp->getSize();
        curvePointSize = (CurvePtNr * (cpSize-1));
        Vector* newCurvePoints = new Vector[curvePointSize];
        for (unsigned i = 0; i < (CurvePtNr * (cpSize - 2)); i++)
        {
            if (curvePoints == NULL) break;
            else newCurvePoints[i] = curvePoints[i];
        }
        delete [] curvePoints;
        curvePoints = newCurvePoints;

        for (unsigned i = 0; i < cp->getSize() -1; i++)
            for (unsigned j = i * CurvePtNr; j < (i+1) * CurvePtNr; j++)
            {
                float ti = cp->getT(i);
                float tip = cp->getT(i + 1);
                float t = ti + (
                        ((tip - ti) / (float)CurvePtNr) * (j - (i * (float)CurvePtNr))
                );
                curvePoints[j] = getCurvePos(i, t);
            }
    }

} dzsCurve(&controlPoints);

class CRSpline: public Curve
{
    static const unsigned pointsBetweenControlPoints = 20;
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
        curvePointSize = points * (cpSize - 1);
        Vector* newCurvePoints = new Vector[curvePointSize];
        delete [] curvePoints;
        curvePoints = newCurvePoints;

        ComputeV();
        for (unsigned i = 0; i < cpSize -1; i++)
            for (unsigned j = i * points; j < (i+1) * points; j++)
            {
                float t = cp->getT(i) + (
                        ((cp->getT(i+1) - cp->getT(i)) / (float)points) * (j - (i * (float)points))
                );
                curvePoints[j] = GetPos(t, i);
            }
    }
} crSpline(&controlPoints);


class BezierCurve : public Curve
{
    static const unsigned resolution = 400;

    // t = [0-1] tartom�ny k�z�tt, i=m=kpMax-1
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
        curvePoints = new Vector[resolution + 1];
        curvePointSize = 0;
    };

    //Bezier nincs sulyozva!
    Vector GetPos(float t)
    {
        return CountBezierPos(t, cp->getSize() - 1, cp->getSize() - 1);
    }

    void computeCurve() {
        float firstT = cp->getT(0);
        float lastT = cp->getT(cp->getSize() - 1);
        float range = lastT - firstT;
        int i = 0;
        for (float t = firstT; t < lastT; t += range / (float) resolution) {
            float t1 = (t - firstT) / range;
            curvePoints[i] = GetPos(t1);
            curvePointSize = i;
            i++;
        }
    }

} bezier(&controlPoints);

class Camera
{
    Vector bottomLeft, topRight;
    Vector originalBottomLeft, originalTopRight;
    bool allowXMove, allowYMove;
public:
    ControlPoint camStart;

    Camera()
    {
        allowXMove = false;
        allowYMove = false;
        originalBottomLeft = bottomLeft = Vector(100.0f, 300.0f);
        originalTopRight = topRight = Vector(200.0f, 400.0f);
    }

    double left(){
        return (double) bottomLeft.x;
    }

    double right(){
        return (double) topRight.x;
    }

    double top(){
        return (double) topRight.y;
    }

    double bottom(){
        return (double) bottomLeft.y;
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

    void startMoving()
    {
        allowXMove = true;
        allowYMove = true;
    }

    void stopXMoving()
    {
        allowXMove = false;
        originalBottomLeft.x = bottomLeft.x;
        originalTopRight.x = topRight.x;
    }

    void stopYMoving()
    {
        allowYMove = false;
        originalBottomLeft.y = bottomLeft.y;
        originalTopRight.y = topRight.y;
    }

    bool isMoving()
    {
        if (allowXMove == true || allowYMove == true)
            return true;
        else return false;
    }

    bool isXMoving()
    {
        return allowXMove;
    }

    bool isYMoving()
    {
        return allowYMove;
    }

} camera;

float ConvertX (int x)
{
    return (float) (camera.left() + ((float)x / screenWidth * 100.0f));
}

float ConvertY (int y)
{
    return (float) (camera.bottom() + (((float)y /screenHeight) * (-1) +1) * 100.0f);
}

// Inicializacio, a program futasanak kezdeten, az OpenGL kontextus letrehozasa utan hivodik meg (ld. main() fv.)
void onInitialization( )
{
    glViewport(0, 0, screenWidth, screenHeight);
    currentSate = ADDING_POINTS;

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
            // calculate curves
        }

    if (currentSate == CAMERA_MOVING)
        for (float ts = tstart; ts < tend; ts += dt) {
            moveCamera(controlPoints.getP(followedControlPoint));
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

    glColor3f(BLACK.r, BLACK.g, BLACK.b);
    for (unsigned i = 0; i < controlPoints.getSize(); i++)
        DrawCircle(controlPoints.getP(i), 2.0f);


    if (controlPoints.getSize() >= 2)
    {
        /*
        glColor3f(BLUE.r, BLUE.g, BLUE.b);
        glBegin(GL_LINE_STRIP);
        for (unsigned i = 0; i < dzsCurve.getCurvePointSize(); i++)
            glVertex2f(dzsCurve.getCurvePoint(i).x, dzsCurve.getCurvePoint(i).y);
        glEnd();
        */

        glColor3f(RED.r, RED.g, RED.b);
        glBegin(GL_LINE_STRIP);
        for (unsigned i = 0; i < bezier.getCurvePointSize(); i++)
            glVertex2f(bezier.getCurvePoint(i).x, bezier.getCurvePoint(i).y);
        glEnd();

        glColor3f(GREEN.r, GREEN.g, GREEN.b);
        glBegin(GL_LINE_STRIP);
        for (unsigned i = 0; i < crSpline.getCurvePointSize(); i++)
            glVertex2f(crSpline.getCurvePoint(i).x, crSpline.getCurvePoint(i).y);
        glEnd();
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
            Vector pos = Vector (ConvertX(x),  ConvertY(y));
            if (controlPoints.getSize() < maxControlPoints)
            {
                controlPoints.add(pos, clickTime/1000.0f);
                if (controlPoints.getSize() >= 2) {
                    dzsCurve.computeCurve();
                    crSpline.computeCurve();
                    bezier.computeCurve();
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
    /*
    if (currentTime - clickTime >= 500 && clickState == MOVING && camera.isMoving() == false)
    {
        Vector CamStartV = Vector(ConvertX(cursorPosX) - ConvertX(clickPosX),
                ConvertY(cursorPosY) - ConvertY(clickPosY)
        );
        CamStartV = CamStartV * -1.0f / 3.0f;
        camera.camStart.v = CamStartV;
        camera.startMoving();
    }
    */

    /*
    if (currentTime - clickTime >= 300)
    {

        switch(clickState)
        {
            case B1UP:
                onB1CLK();
                clickState = START;
                break;
            case B2UP:
                onB2CLK();
                clickState = START;
                break;
            default: break;
        }
    }
    */

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
