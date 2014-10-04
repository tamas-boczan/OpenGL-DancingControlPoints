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

// színek
const Color BLACK = Color(0.0f, 0.0f, 0.0f);
const Color GREY = Color(0.85f, 0.85f, 0.85f);
// Türkíz szín forrása: www.opengl.org/discussion_boards/showthread.php/132502-Color-tables
const Color TURQUOISE = Color(0.678f, 0.918f, 0.918f);
const Color RED = Color(1.0f, 0.0f, 0.0f);
const Color GREEN = Color(0.0f, 0.6f, 0.1f);


long currentTime;
long clickTime;
long circulationStartTime;
unsigned int followedControlPoint;

enum possibleStates{ADDING_POINTS, CIRCULATE, CAMERA_MOVING} currentSate;

struct ControlPoint
{
    Vector p;
    Vector a;
    Vector v;
    float t;
};

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

class Curve
{
public:
    static const unsigned MaxCtrlPts = 100;
    static const unsigned CurvePtNr = 20;

private:

    ControlPoint* cp;
    Vector* CurvePoints;
    size_t cpSize;

    Vector a0(int i)
    {
        return cp[i].p;
    }

    Vector a1(int i)
    {
        return cp[i].v;
    }

    Vector a2(int i)
    {
        return cp[i].a / 2.0f;
    }

    Vector a3(int i)
    {

        return (
                ((cp[i+1].p * 4.0f - cp[i].p * 4.0f) / (pow((cp[i+1].t - cp[i].t), 3)))
                        - ((cp[i+1].v + cp[i].v * 3.0f) / (pow((cp[i+1].t - cp[i].t), 2)))
                        - (cp[i].a / (cp[i+1].t - cp[i].t))
        );
    }

    Vector a4(int i)
    {
        return (
                (((cp[i+1].p * (-3.0f)) + (cp[i].p * 3.0f)) / (pow((cp[i+1].t - cp[i].t), 4)))
                        + ((cp[i+1].v + (cp[i].v * 2.0f)) / (pow((cp[i+1].t - cp[i].t),3)))
                        + ((cp[i].a *0.5f) / (pow((cp[i+1].t - cp[i].t), 2)))
        );
    }

    void computeV()
    {
        cp[cpSize-1].v = Vector(0.0f, 0.0f, 0.0f);
        for (unsigned i = 1; i < cpSize - 1; i++)
            cp[i].v = (
                    ((cp[i].p - cp[i-1].p) / (cp[i].t - cp[i-1].t)) +
                            ((cp[i+1].p - cp[i].p) / (cp[i+1].t - cp[i].t))
            ) / 2.0f;
    }

    void computeA()
    {
        for (unsigned i = 0; i < cpSize -1; i++)
        {

            cp[i+1].a = (a4(i) * 12.0f * (pow((cp[i+1].t - cp[i].t),2))) +
                    (a3(i) *  6.0f * (cp[i+1].t - cp[i].t)) +
                    (a2(i) *  2.0f);
        }
    }


    int findCP(float t)
    {
        for(unsigned i = 1; i < cpSize; i++)
            if (cp[i].t > t)
                return i-1;
        return cpSize-1;
    }
public:
    Curve()
    {
        cp = NULL;
        CurvePoints = NULL;
        cpSize=0;
    }

    ~Curve()
    {
        delete[] cp;
        delete[] CurvePoints;
    }

    size_t getCpNr()
    {
        return cpSize;
    }

    Vector getP(unsigned i)
    {
        return  cp[i].p;
    }

    float getT(unsigned i)
    {
        return cp[i].t;
    }

    float getLastT()
    {
        return cp[cpSize-1].t;
    }

    void setP(unsigned i, Vector Position)
    {
        cp[i].p = Position;
    }

    void setT(unsigned i, float weight)
    {
        cp[i].t = weight;
    }

    void setV(unsigned i, Vector Velocity)
    {
        cp[i].v = Velocity;
    }

    void setA(unsigned i, Vector Acceleration)
    {
        cp[i].a = Acceleration;
    }

    void addCP(Vector Position, float weight)
    {
        ControlPoint* newcp = new ControlPoint[cpSize+1];
        for (unsigned i = 0; i < cpSize; i++)
            newcp[i] = cp[i];
        newcp[cpSize].p = Position;
        newcp[cpSize].t = weight;
        delete [] cp;
        cp = newcp;

        cpSize++;
        if (cpSize >= 3)
        {
            Vector* newCurvePoints = new Vector[(CurvePtNr * (cpSize-1))];
            for (unsigned i = 0; i < (CurvePtNr * (cpSize - 2)); i++)
            {
                if (CurvePoints == NULL) break;
                else newCurvePoints[i] = CurvePoints[i];
            }
            delete [] CurvePoints;
            CurvePoints = newCurvePoints;
        }
    }

    Vector getCurvePos(int previousControlPointIndex, float weight)
    {
        int i = previousControlPointIndex;
        float t = weight;

        return	(a4(i) * pow((t-cp[i].t),4)) +
                (a3(i) * pow((t-cp[i].t),3)) +
                (a2(i) * pow((t-cp[i].t),2)) +
                (a1(i) * (t-cp[i].t)) +
                a0(i);
    }

    Vector getCurvePos(float weight)
    {
        int previousControlPointIndex = findCP(weight);
        return getCurvePos(previousControlPointIndex, weight);
    }

    Vector getCurvePoint(int i)
    {
        return CurvePoints[i];
    }

    void computeCurve()
    {
        computeV();
        computeA();


        for (unsigned i = 0; i < cpSize -1; i++)
            for (unsigned j = i * CurvePtNr; j < (i+1) * CurvePtNr; j++)
            {
                float t = cp[i].t + (
                        ((cp[i+1].t - cp[i].t) / (float)CurvePtNr) * (j - (i * (float)CurvePtNr))
                );
                CurvePoints[j] = getCurvePos(i, t);

            }

    }
} DzsCurve;

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
            moveCamera(DzsCurve.getP(followedControlPoint));
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
    for (unsigned i = 0; i < DzsCurve.getCpNr(); i++)
        DrawCircle(DzsCurve.getP(i), 2.0f);

    glColor3f(GREEN.r, GREEN.g, GREEN.b);
    if (DzsCurve.getCpNr() >=3)
    {
        glBegin(GL_LINE_STRIP);
        for (unsigned i = 0; i < (DzsCurve.getCpNr()-1) * DzsCurve.CurvePtNr; i++)
            glVertex2f(DzsCurve.getCurvePoint(i).x, DzsCurve.getCurvePoint(i).y);
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
            if (DzsCurve.getCpNr() < DzsCurve.MaxCtrlPts)
            {
                DzsCurve.addCP(pos, clickTime/1000.0f);
                if (DzsCurve.getCpNr() >= 3)
                    DzsCurve.computeCurve();
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
