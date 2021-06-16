#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <stdlib.h>
#include <stdio.h>
#include <windows.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include<bits/stdc++.h>
using namespace std;

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

using namespace std;

#define PI 3.14159265

const int width = 800;
const int height = 500;
const float rat = 1.0*width/height;
float spt_cutoff = 30;
int anglex= 0, angley = 0, anglez = 0;
int window;
int wired=0;
int shcpt=1;
int animat = 0;
const int L=20;
const int dgre=3;
int ncpt=L+1;
int clikd=0;
const int nt = 40;				//number of slices along x-direction
const int ntheta = 20;
bool flagg=false,flagg2=false,flagg3=false,flagg5=false,flagg6=false,flagg7=false,flagg8=false,flagg9=false,flagg10=false;

float spot_cutoff = 40;
bool spot_light= false;

GLfloat eyeX = 0;
GLfloat eyeY = 6;
GLfloat eyeZ = 38;

GLfloat lookX = 0;
GLfloat lookY = 6;
GLfloat lookZ = 0;

float rot = 0;
float fan_rot = 0;
int wndw_rot=0;
float inter_wndw_rot=0;


GLboolean l_zero=false,l_one=false;

bool target_on=true,fan_on=true;


///light 0
GLfloat l_no[] = {0, 0, 0, 1.0};
GLfloat l_amb[] = {0.5, 0.5, 0.5, 1.0};
GLfloat l_dif[] = {1,1,1,1};
GLfloat l_spec[] = {.6,.6,.6,1};
GLfloat l_pos[] = {6.5,10,-4,1.0};


///light 1
GLfloat l_no1[] = {0, 0, 0, 1.0};
GLfloat l_amb1[] = {0.3, 0.3, 0.3, 1.0};
GLfloat l_dif1[] = {1,1,1,1};
GLfloat l_spec1[] = {.5,.5,.5,1};
GLfloat l_pos1[] = {-9.5,10,-4,1.0};

///light 2
GLfloat l_no2[] = {0, 0, 0, 1.0};
GLfloat l_amb2[] = {0.3, 0.3, 0.3, 1};
GLfloat l_dif2[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat l_spec2[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat l_pos2[] = {1.5,2.5,18,1.0 };

float theata=0;
float z=2.8,x=4,c=1;
bool flagg4=false;

unsigned int ID,ID2;

GLfloat ctrlpoints[L+1][3] =
{
    { 0.0, 0.0, 0.0},{0.1,1.4,0},{1,2.4,0},{0.8,3.0,0},{1.3,4,0},{1.8,4.5,0}

};

GLfloat ctrlpoints1[L+1][3] =
{
     { 0.0, 0.0, 0.0}, { -0.3, 0.5, 0.0},
    { 0.1, 1.7, 0.0},{ 0.5, 1.5, 0.0},
    {1.0, 1.5, 0.0}, {1.4, 1.4, 0.0},
    {1.8, 0.4, 0.0},{2.2, 0.4, 0.0},
    {2.6, 1.5, 0.0}, {3.0, 1.4, 0.0},
    {3.4, 1.4, 0.0},{3.8, 1.4, 0.0},
    {4.2, 1.0, 0.0},{4.6, 1.0, 0.0},
    {5.0, 1.0, 0.0},{5.4, 1.0, 0.0},
    {5.8, 0.5, 0.0},{6.2, 0.5, 0.0},
    {6.6, 0.5, 0.0},{7.2, 0.2, 0.0},
    {6.8, 0.52, 0.0}

};
float wcsClkDn[3],wcsClkUp[3];

class point1
{
public:
    point1()
    {
        x=0;
        y=0;
    }
    int x;
    int y;
} clkpt[2];
int flag=0;
GLint viewport[4]; //var to hold the viewport info
GLdouble modelview[16]; //var to hold the modelview info
GLdouble projection[16];

class BmpLoader
{
public:
    unsigned char* textureData;
    int iWidth, iHeight;

    BmpLoader(const char* filename)
    {
        FILE *file=0;
        file=fopen(filename, "rb");
        if(!file)
            std::cout<<"File not found"<<std::endl;
        fread(&bfh, sizeof(BITMAPFILEHEADER),1,file);
        if(bfh.bfType != 0x4D42)
            std::cout<<"Not a valid bitmap"<<std::endl;
        fread(&bih, sizeof(BITMAPINFOHEADER),1,file);
        if(bih.biSizeImage==0)
            bih.biSizeImage=bih.biHeight*bih.biWidth*3;
        textureData = new unsigned char[bih.biSizeImage];
        fseek(file, bfh.bfOffBits, SEEK_SET);
        fread(textureData, 1, bih.biSizeImage, file);
        unsigned char temp;
        for(int i=0; i<bih.biSizeImage; i+=3)
        {
            temp = textureData[i];
            textureData[i] = textureData[i+2];
            textureData[i+2] = temp;

        }

        iWidth = bih.biWidth;
        iHeight = bih.biHeight;
        fclose(file);
    }
    ~BmpLoader()
    {
        delete [] textureData;
    }

private:
    BITMAPFILEHEADER bfh;
    BITMAPINFOHEADER bih;
};
void scsToWcs(float sx,float sy, float wcsv[3] );
void processMouse(int button, int state, int x, int y);

///////////////////////////
void scsToWcs(float sx,float sy, float wcsv[3] )
{

    GLfloat winX, winY, winZ; //variables to hold screen x,y,z coordinates
    GLdouble worldX, worldY, worldZ; //variables to hold world x,y,z coordinates

    //glGetDoublev( GL_MODELVIEW_MATRIX, modelview ); //get the modelview info
    glGetDoublev( GL_PROJECTION_MATRIX, projection ); //get the projection matrix info
    glGetIntegerv( GL_VIEWPORT, viewport ); //get the viewport info

    winX = sx;
    winY = (float)viewport[3] - (float)sy;
    winZ = 0;

    //get the world coordinates from the screen coordinates
    gluUnProject( winX, winY, winZ, modelview, projection, viewport, &worldX, &worldY, &worldZ);
    wcsv[0]=worldX;
    wcsv[1]=worldY;
    wcsv[2]=worldZ;


}
void processMouse(int button, int state, int x, int y)
{
    if(button==GLUT_LEFT_BUTTON && state==GLUT_DOWN)
    {
        if(flag!=1)
        {
            flag=1;
            clkpt[0].x=x;
            clkpt[0].y=y;
        }


        scsToWcs(clkpt[0].x,clkpt[0].y,wcsClkDn);
        cout<<"\nD: "<<x<<" "<<y<<" wcs: "<<wcsClkDn[0]<<" "<<wcsClkDn[1];
    }
    else if(button==GLUT_LEFT_BUTTON && state==GLUT_UP)
    {
        if (flag==1)
        {
            clkpt[1].x=x;
            clkpt[1].y=y;
            flag=0;
        }
        float wcs[3];
        scsToWcs(clkpt[1].x,clkpt[1].y,wcsClkUp);
        cout<<"\nU: "<<x<<" "<<y<<" wcs: "<<wcsClkUp[0]<<" "<<wcsClkUp[1];

        clikd=!clikd;
    }
}

long long nCr(int n, int r)
{
    if(r > n / 2) r = n - r; // because C(n, r) == C(n, n - r)
    long long ans = 1;
    int i;

    for(i = 1; i <= r; i++)
    {
        ans *= n - r + i;
        ans /= i;
    }

    return ans;
}


void BezierCurve ( double t,  float xy[2])
{
    double y=0;
    double x=0;
    t=t>1.0?1.0:t;
    for(int i=0; i<=L; i++)
    {
        int ncr=nCr(L,i);
        double oneMinusTpow=pow(1-t,double(L-i));
        double tPow=pow(t,double(i));
        double coef=oneMinusTpow*tPow*ncr;
        x+=coef*ctrlpoints[i][0];
        y+=coef*ctrlpoints[i][1];

    }
    xy[0] = float(x);
    xy[1] = float(y);

    //return y;
}

void BezierCurve1 ( double t,  float xy[2])
{
    double y=0;
    double x=0;
    t=t>1.0?1.0:t;
    for(int i=0; i<=L; i++)
    {
        int ncr=nCr(L,i);
        double oneMinusTpow=pow(1-t,double(L-i));
        double tPow=pow(t,double(i));
        double coef=oneMinusTpow*tPow*ncr;
        x+=coef*ctrlpoints1[i][0];
        y+=coef*ctrlpoints1[i][1];

    }
    xy[0] = float(x);
    xy[1] = float(y);

    //return y;
}

void setNormal(GLfloat x1, GLfloat y1,GLfloat z1, GLfloat x2, GLfloat y2,GLfloat z2, GLfloat x3, GLfloat y3,GLfloat z3)
{
    GLfloat Ux, Uy, Uz, Vx, Vy, Vz, Nx, Ny, Nz;

    Ux = x2-x1;
    Uy = y2-y1;
    Uz = z2-z1;

    Vx = x3-x1;
    Vy = y3-y1;
    Vz = z3-z1;

    Nx = Uy*Vz - Uz*Vy;
    Ny = Uz*Vx - Ux*Vz;
    Nz = Ux*Vy - Uy*Vx;

    glNormal3f(-Nx,-Ny,-Nz);
}

void bottleBezier()
{
    int i, j;
    float x, y, z, r;				//current coordinates
    float x1, y1, z1, r1;			//next coordinates
    float theta;

    const float startx = 0, endx = ctrlpoints[L][0];
    //number of angular slices
    const float dx = (endx - startx) / nt;	//x step size
    const float dtheta = 2*PI / ntheta;		//angular step size

    float t=0;
    float dt=1.0/nt;
    float xy[2];
    BezierCurve( t,  xy);
    x = xy[0];
    r = xy[1];
    //rotate about z-axis
    float p1x,p1y,p1z,p2x,p2y,p2z;
    for ( i = 0; i < nt; ++i )  			//step through x
    {
        theta = 0;
        t+=dt;
        BezierCurve( t,  xy);
        x1 = xy[0];
        r1 = xy[1];

        //draw the surface composed of quadrilaterals by sweeping theta
        glBegin( GL_QUAD_STRIP );
        ///glBegin( GL_QUADS );
        for ( j = 0; j <= ntheta; ++j )
        {
            theta += dtheta;
            double cosa = cos( theta );
            double sina = sin ( theta );
            y = r * cosa;
            y1 = r1 * cosa;	//current and next y
            z = r * sina;
            z1 = r1 * sina;	//current and next z

            //edge from point at x to point at next x
            glVertex3f (x, y, z);

            if(j>0)
            {
                setNormal(p1x,p1y,p1z,p2x,p2y,p2z,x, y, z);
            }
            else
            {
                p1x=x;
                p1y=y;
                p1z=z;
                p2x=x1;
                p2y=y1;
                p2z=z1;

            }
            glVertex3f (x1, y1, z1);

            //forms quad with next pair of points with incremented theta value
        }
        glEnd();
        x = x1;
        r = r1;
    } //for i

}


///bottle
void bottleBezier1()
{
    int i, j;
    float x, y, z, r;				//current coordinates
    float x1, y1, z1, r1;			//next coordinates
    float theta;

    const float startx = 0, endx = ctrlpoints1[L][0];
    //number of angular slices
    const float dx = (endx - startx) / nt;	//x step size
    const float dtheta = 2*PI / ntheta;		//angular step size

    float t=0;
    float dt=1.0/nt;
    float xy[2];
    BezierCurve1( t,  xy);
    x = xy[0];
    r = xy[1];
    //rotate about z-axis
    float p1x,p1y,p1z,p2x,p2y,p2z;
    for ( i = 0; i < nt; ++i )  			//step through x
    {
        theta = 0;
        t+=dt;
        BezierCurve1( t,  xy);
        x1 = xy[0];
        r1 = xy[1];

        //draw the surface composed of quadrilaterals by sweeping theta
        glBegin( GL_QUAD_STRIP );
        ///glBegin( GL_QUADS );
        for ( j = 0; j <= ntheta; ++j )
        {
            theta += dtheta;
            double cosa = cos( theta );
            double sina = sin ( theta );
            y = r * cosa;
            y1 = r1 * cosa;	//current and next y
            z = r * sina;
            z1 = r1 * sina;	//current and next z

            //edge from point at x to point at next x
            glVertex3f (x, y, z);

            if(j>0)
            {
                setNormal(p1x,p1y,p1z,p2x,p2y,p2z,x, y, z);
            }
            else
            {
                p1x=x;
                p1y=y;
                p1z=z;
                p2x=x1;
                p2y=y1;
                p2z=z1;

            }
            glVertex3f (x1, y1, z1);

            //forms quad with next pair of points with incremented theta value
        }
        glEnd();
        x = x1;
        r = r1;
    } //for i

}

void showControlPoints()
{
    glPointSize(5.0);
    glColor3f(1.0, 0.0, 1.0);
    glBegin(GL_POINTS);
    for (int i = 0; i <=L; i++)
        glVertex3fv(&ctrlpoints[i][0]);
    glEnd();
}


static GLfloat v_cube[8][3] =
{
    {0,0,0},
    {0,0,1},
    {0,1,0},
    {0,1,1},

    {1,0,0},
    {1,0,1},
    {1,1,0},
    {1,1,1}
};

static GLubyte c_ind[6][4] =
{
    {3,1,5,7},  //front
    {6,4,0,2},  //back
    {2,3,7,6},  //top
    {1,0,4,5},  //bottom
    {7,5,4,6},  //right
    {2,0,1,3}   //left
};

static void getNormal3p(GLfloat x1, GLfloat y1, GLfloat z1,
                        GLfloat x2, GLfloat y2, GLfloat z2,
                        GLfloat x3, GLfloat y3, GLfloat z3)
{
    GLfloat Ux, Uy, Uz, Vx, Vy, Vz, Nx, Ny, Nz;

    Ux = x2-x1;
    Uy = y2-y1;
    Uz = z2-z1;

    Vx = x3-x1;
    Vy = y3-y1;
    Vz = z3-z1;

    Nx = Uy*Vz - Uz*Vy;
    Ny = Uz*Vx - Ux*Vz;
    Nz = Ux*Vy - Uy*Vx;

    glNormal3f(Nx,Ny,Nz);
}

void set_mat_prop(float colR=.3, float colG=.3, float colB=.3, bool em=false, float shine=60)
{
    GLfloat no_mat[] = { 0.0, 0.0, 0.0, 1.0 };
    GLfloat mat_ambient[] = { colR, colG, colB, 1.0 };
    GLfloat mat_diffuse[] = { colR, colG, colB, 1.0 };
    GLfloat mat_specular[] = { 0.5, 0.5, 0.5, 1.0 };
    GLfloat mat_emission[] = {colR, colG, colB, 1.0};
    GLfloat mat_shininess[] = {shine};

    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient);
    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse);
    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess);

    if(em)
        glMaterialfv( GL_FRONT, GL_EMISSION, mat_emission);
    else
        glMaterialfv( GL_FRONT, GL_EMISSION, no_mat);
}

void cube(float colR=.5, float colG=.5, float colB=.5,
          bool em=true, float shine=60)
{
    set_mat_prop(colR,colG,colB,em,shine);

    glBegin(GL_QUADS);
    for (GLint i = 0; i <6; i++)
    {
        getNormal3p(v_cube[c_ind[i][0]][0], v_cube[c_ind[i][0]][1], v_cube[c_ind[i][0]][2],
                    v_cube[c_ind[i][1]][0], v_cube[c_ind[i][1]][1], v_cube[c_ind[i][1]][2],
                    v_cube[c_ind[i][2]][0], v_cube[c_ind[i][2]][1], v_cube[c_ind[i][2]][2]);

        glTexCoord2f(0,1);
        glVertex3fv(&v_cube[c_ind[i][0]][0]);
        glTexCoord2f(0,0);
        glVertex3fv(&v_cube[c_ind[i][1]][0]);
        glTexCoord2f(1,0);
        glVertex3fv(&v_cube[c_ind[i][2]][0]);
        glTexCoord2f(1,1);
        glVertex3fv(&v_cube[c_ind[i][3]][0]);
    }
    glEnd();

}




GLubyte* make_bw_tiles_texture(int tile_width, int tile_height, int tex_width=2048, int tex_height=2048)
{
    GLubyte* bw_tiles = new GLubyte[tex_width*tex_height*3];

    for(int i=0; i<tex_height; i++)
    {
        for(int j=0; j<tex_width; j++)
        {
            int c=(i/tile_width + j/tile_width) % 2;

            if(c==0)
            {
                bw_tiles[(i*tex_width+j)*3+0]=255;
                bw_tiles[(i*tex_width+j)*3+1]=255;
                bw_tiles[(i*tex_width+j)*3+2]=255;
            }
            else
            {
                bw_tiles[(i*tex_width+j)*3+0]=0;
                bw_tiles[(i*tex_width+j)*3+1]=0;
                bw_tiles[(i*tex_width+j)*3+2]=0;
            }
        }
    }

    return bw_tiles;
}



void LoadTexture(const char*filename)
{
    glGenTextures(1, &ID);
    glBindTexture(GL_TEXTURE_2D, ID);
    glPixelStorei(GL_UNPACK_ALIGNMENT, ID);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT );
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    BmpLoader bl(filename);
    gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGB, bl.iWidth, bl.iHeight, GL_RGB, GL_UNSIGNED_BYTE, bl.textureData );


}

void animate()
{
    if (flagg == true)
    {
        theata-= 0.1;
        if(theata < -180)
        {
            theata = -180;
            flagg=false;
        }
    }
    glutPostRedisplay();
    if(flagg2==true)
    {
        theata+=0.1;
        if(theata>0)
        {
            theata=0;
            flagg2=false;
        }
    }
    glutPostRedisplay();

    if(flagg3==true)
    {
        z-=0.1;
        if(z<-1)
        {
            z=-1;
            flagg3=false;
        }
    }
    glutPostRedisplay();

    if(flagg4==true)
    {
        z+=0.1;
        if(z>2.8)
        {
            z=2.8;
            flagg4=false;
        }
    }
    glutPostRedisplay();

    if(flagg5==true)
    {
        x-=0.1;
        if(x<-4)
        {
            flagg5=false;
            x=-4;
        }
    }
    glutPostRedisplay();

    if(flagg6==true)
    {
        x+=0.1;
        if(x>4)
        {
            flagg5=false;
            x=4;
        }
    }
    glutPostRedisplay();

    if(flagg7==true)
    {

        x+=0.1;
        z+=0.047;
        if(x>4 && z>2.8)
        {
            flagg7=false;

            x=4;
            z=2.8;
        }

    }
    glutPostRedisplay();

    if(flagg8==true)
    {

        x-=0.01;
        z-=0.0047;
        if(x>-4 && z>-1)
        {
            flagg8=false;

            x=-4;
            z=-1;
        }

    }
    glutPostRedisplay();

    if(flagg9)
    {
        c-=0.1;
        if(c<-12)
        {
            flagg9=false;
            c=-12;
        }
    }
    glutPostRedisplay();

    if(flagg10)
    {
        c+=0.1;
        if(c>1)
        {
            flagg10=false;
            c=1;
        }
    }
    glutPostRedisplay();
}
void ball()
{

    glPushMatrix();
    set_mat_prop(0,0,0,true,60);

    glTranslated(x,.3,z);

    glutSolidSphere(0.2,20,20);
    glPopMatrix();
}
void table_chair()
{

}
void wall()
{
    ///left wall
    glPushMatrix();
    glTranslated(-12,0,-12);
    glScaled(.5,12,19);
    cube(0.627, 0.322, 0.176);
    glPopMatrix();

    ///right wall
    glPushMatrix();
    glTranslated(12,0,-12);
    glScaled(.5,12,19);
    cube(0.627, 0.322, 0.176);
    glPopMatrix();

    /// back wall

    glPushMatrix();
    glTranslated(-12,0,-12);
    glScaled(24,12,.5);
    cube(0.627, 0.322, 0.176);
    glPopMatrix();

    ///front wall

    glPushMatrix();
    glTranslated(-12,0,7);
    glScaled(9,12,.5);
    cube(0.627, 0.322, 0.176);
    glPopMatrix();

    glPushMatrix();
    glTranslated(5.5,0,7);
    glScaled(9,12,.5);
    cube(0.627, 0.322, 0.176);
    glPopMatrix();

    glPushMatrix();
    glTranslated(-3,0,7);
    glRotated(theata,0,1,0);


    glScaled(9,12,.2);

    cube(0.561, 0.737, 0.561);
    glPopMatrix();




}
void build()
{
    glPushMatrix();
    glTranslated(-12,0,7.5);
    glScaled(9,.1,17);
    cube(1,1,1);
    glPopMatrix();


    glPushMatrix();
    glTranslated(5.5,0,7.5);
    glScaled(9,.1,17);
    cube(1.000, 0.894, 0.769);
    glPopMatrix();

    ///another building1
    glPushMatrix();
    glScaled(.5,.5,.5);
    glTranslated(-32,0,20);

    glPushMatrix();
    glTranslated(8,0,5);
    glScaled(2,16,1.5);
    cube(1.000, 0.627, 0.478);
    glPopMatrix();

    glPushMatrix();
    glTranslated(8,0,8);
    glScaled(2,16,3);
    cube(1.000, 0.627, 0.478);
    glPopMatrix();

    glPushMatrix();
    glTranslated(8,0,12);
    glScaled(2,16,1.5);
    cube(1.000, 0.627, 0.478);
    glPopMatrix();

    glPushMatrix();
    glTranslated(8,16,5);
    glScaled(2,2,8.5);
    cube(0.184, 0.310, 0.310);
    glPopMatrix();

    glPushMatrix();
    glTranslated(8,10,6.5);
    glScaled(2,2,1.5);
    cube(0.184, 0.310, 0.310);
    glPopMatrix();

    glPushMatrix();
    glTranslated(8,10,11);
    glScaled(2,2,1);
    cube(0.184, 0.310, 0.310);
    glPopMatrix();

    glPushMatrix();
    glTranslated(9,12,11);
    glScaled(.5,4,1);
    cube(0.941, 0.902, 0.549);
    glPopMatrix();

    glPushMatrix();
    glTranslated(9,12,6.5);
    glScaled(.5,4,1.5);
    cube(0.941, 0.902, 0.549);
    glPopMatrix();


    glPushMatrix();
    glTranslated(8,4,11);
    glScaled(2,2,1);
    cube(0.184, 0.310, 0.310);
    glPopMatrix();

    glPushMatrix();
    glTranslated(8,4,6.5);
    glScaled(2,2,1.5);
    cube(0.184, 0.310, 0.310);
    glPopMatrix();

    glPushMatrix();
    glTranslated(9,6,11);
    glScaled(.5,4,1);
    cube(0.941, 0.902, 0.549);
    glPopMatrix();

    glPushMatrix();
    glTranslated(9,6,6.5);
    glScaled(.5,4,1.5);
    cube(0.941, 0.902, 0.549);
    glPopMatrix();

    glPushMatrix();
    glTranslated(9,0,11);
    glScaled(.5,4,1);
    cube(0.941, 0.902, 0.549);
    glPopMatrix();

    glPushMatrix();
    glTranslated(9,0,6.5);
    glScaled(.5,4,1.5);
    cube(0.941, 0.902, 0.549);
    glPopMatrix();

    glPushMatrix();
    glTranslated(8,0,6.5);
    glScaled(.5,.5,1.5);
    cube(0.941, 0.902, 0.549);
    glPopMatrix();

    glPushMatrix();
    glTranslated(8,0,11);
    glScaled(.5,.5,1.5);
    cube(0.941, 0.902, 0.549);
    glPopMatrix();

    glPopMatrix();


    ///build 2

    glPushMatrix();
    glScaled(.5,.5,.5);
    glTranslated(-32,0,32);

    glPushMatrix();
    glTranslated(8,0,5);
    glScaled(2,16,1.5);
    cube(0.941, 0.502, 0.502);
    glPopMatrix();

    glPushMatrix();
    glTranslated(8,0,8);
    glScaled(2,16,3);
    cube(0.941, 0.502, 0.502);
    glPopMatrix();

    glPushMatrix();
    glTranslated(8,0,12);
    glScaled(2,16,1.5);
    cube(0.941, 0.502, 0.502);
    glPopMatrix();

    glPushMatrix();
    glTranslated(8,16,5);
    glScaled(2,2,8.5);
    cube(0.184, 0.310, 0.310);
    glPopMatrix();

    glPushMatrix();
    glTranslated(8,10,6.5);
    glScaled(2,2,1.5);
    cube(0.184, 0.310, 0.310);
    glPopMatrix();

    glPushMatrix();
    glTranslated(8,10,11);
    glScaled(2,2,1);
    cube(0.184, 0.310, 0.310);
    glPopMatrix();

    glPushMatrix();
    glTranslated(9,12,11);
    glScaled(.5,4,1);
    cube(0.941, 0.902, 0.549);
    glPopMatrix();

    glPushMatrix();
    glTranslated(9,12,6.5);
    glScaled(.5,4,1.5);
    cube(0.941, 0.902, 0.549);
    glPopMatrix();


    glPushMatrix();
    glTranslated(8,4,11);
    glScaled(2,2,1);
    cube(0.184, 0.310, 0.310);
    glPopMatrix();

    glPushMatrix();
    glTranslated(8,4,6.5);
    glScaled(2,2,1.5);
    cube(0.184, 0.310, 0.310);
    glPopMatrix();

    glPushMatrix();
    glTranslated(9,6,11);
    glScaled(.5,4,1);
    cube(0.941, 0.902, 0.549);
    glPopMatrix();

    glPushMatrix();
    glTranslated(9,6,6.5);
    glScaled(.5,4,1.5);
    cube(0.941, 0.902, 0.549);
    glPopMatrix();

    glPushMatrix();
    glTranslated(9,0,11);
    glScaled(.5,4,1);
    cube(0.941, 0.902, 0.549);
    glPopMatrix();

    glPushMatrix();
    glTranslated(9,0,6.5);
    glScaled(.5,4,1.5);
    cube(0.941, 0.902, 0.549);
    glPopMatrix();

    glPushMatrix();
    glTranslated(8,0,6.5);
    glScaled(.5,.5,1.5);
    cube(0.941, 0.902, 0.549);
    glPopMatrix();

    glPushMatrix();
    glTranslated(8,0,11);
    glScaled(.5,.5,1.5);
    cube(0.941, 0.902, 0.549);
    glPopMatrix();

    glPopMatrix();


    ///build 3

    glPushMatrix();
    glScaled(.5,.5,.5);
    glTranslated(18,0,32);

    glPushMatrix();
    glTranslated(8,0,5);
    glScaled(2,16,1.5);
    cube(0.647, 0.165, 0.165);
    glPopMatrix();

    glPushMatrix();
    glTranslated(8,0,8);
    glScaled(2,16,3);
    cube(0.647, 0.165, 0.165);
    glPopMatrix();

    glPushMatrix();
    glTranslated(8,0,12);
    glScaled(2,16,1.5);
    cube(0.647, 0.165, 0.165);
    glPopMatrix();

    glPushMatrix();
    glTranslated(8,16,5);
    glScaled(2,2,8.5);
    cube(0.184, 0.310, 0.310);
    glPopMatrix();

    glPushMatrix();
    glTranslated(8,10,6.5);
    glScaled(2,2,1.5);
    cube(0.184, 0.310, 0.310);
    glPopMatrix();

    glPushMatrix();
    glTranslated(8,10,11);
    glScaled(2,2,1);
    cube(0.184, 0.310, 0.310);
    glPopMatrix();

    glPushMatrix();
    glTranslated(9,12,11);
    glScaled(.5,4,1);
    cube(0.941, 0.902, 0.549);
    glPopMatrix();

    glPushMatrix();
    glTranslated(9,12,6.5);
    glScaled(.5,4,1.5);
    cube(0.941, 0.902, 0.549);
    glPopMatrix();


    glPushMatrix();
    glTranslated(8,4,11);
    glScaled(2,2,1);
    cube(0.184, 0.310, 0.310);
    glPopMatrix();

    glPushMatrix();
    glTranslated(8,4,6.5);
    glScaled(2,2,1.5);
    cube(0.184, 0.310, 0.310);
    glPopMatrix();

    glPushMatrix();
    glTranslated(9,6,11);
    glScaled(.5,4,1);
    cube(0.941, 0.902, 0.549);
    glPopMatrix();

    glPushMatrix();
    glTranslated(9,6,6.5);
    glScaled(.5,4,1.5);
    cube(0.941, 0.902, 0.549);
    glPopMatrix();

    glPushMatrix();
    glTranslated(9,0,11);
    glScaled(.5,4,1);
    cube(0.941, 0.902, 0.549);
    glPopMatrix();

    glPushMatrix();
    glTranslated(9,0,6.5);
    glScaled(.5,4,1.5);
    cube(0.941, 0.902, 0.549);
    glPopMatrix();

    glPushMatrix();
    glTranslated(8,0,6.5);
    glScaled(.5,.5,1.5);
    cube(0.941, 0.902, 0.549);
    glPopMatrix();

    glPushMatrix();
    glTranslated(8,0,11);
    glScaled(.5,.5,1.5);
    cube(0.941, 0.902, 0.549);
    glPopMatrix();

    glPopMatrix();


    ///build 4

    glPushMatrix();
    glScaled(.5,.5,.5);
    glTranslated(18,0,20);

    glPushMatrix();
    glTranslated(8,0,5);
    glScaled(2,16,1.5);
    cube(0.941, 1.000, 0.941);
    glPopMatrix();

    glPushMatrix();
    glTranslated(8,0,8);
    glScaled(2,16,3);
    cube(0.941, 1.000, 0.941);
    glPopMatrix();

    glPushMatrix();
    glTranslated(8,0,12);
    glScaled(2,16,1.5);
    cube(0.941, 1.000, 0.941);
    glPopMatrix();

    glPushMatrix();
    glTranslated(8,16,5);
    glScaled(2,2,8.5);
    cube(0.184, 0.310, 0.310);
    glPopMatrix();

    glPushMatrix();
    glTranslated(8,10,6.5);
    glScaled(2,2,1.5);
    cube(0.184, 0.310, 0.310);
    glPopMatrix();

    glPushMatrix();
    glTranslated(8,10,11);
    glScaled(2,2,1);
    cube(0.184, 0.310, 0.310);
    glPopMatrix();

    glPushMatrix();
    glTranslated(9,12,11);
    glScaled(.5,4,1);
    cube(0.941, 0.902, 0.549);
    glPopMatrix();

    glPushMatrix();
    glTranslated(9,12,6.5);
    glScaled(.5,4,1.5);
    cube(0.941, 0.902, 0.549);
    glPopMatrix();


    glPushMatrix();
    glTranslated(8,4,11);
    glScaled(2,2,1);
    cube(0.184, 0.310, 0.310);
    glPopMatrix();

    glPushMatrix();
    glTranslated(8,4,6.5);
    glScaled(2,2,1.5);
    cube(0.184, 0.310, 0.310);
    glPopMatrix();

    glPushMatrix();
    glTranslated(9,6,11);
    glScaled(.5,4,1);
    cube(0.941, 0.902, 0.549);
    glPopMatrix();

    glPushMatrix();
    glTranslated(9,6,6.5);
    glScaled(.5,4,1.5);
    cube(0.941, 0.902, 0.549);
    glPopMatrix();

    glPushMatrix();
    glTranslated(9,0,11);
    glScaled(.5,4,1);
    cube(0.941, 0.902, 0.549);
    glPopMatrix();

    glPushMatrix();
    glTranslated(9,0,6.5);
    glScaled(.5,4,1.5);
    cube(0.941, 0.902, 0.549);
    glPopMatrix();

    glPushMatrix();
    glTranslated(8,0,6.5);
    glScaled(.5,.5,1.5);
    cube(0.941, 0.902, 0.549);
    glPopMatrix();

    glPushMatrix();
    glTranslated(8,0,11);
    glScaled(.5,.5,1.5);
    cube(0.941, 0.902, 0.549);
    glPopMatrix();

    glPopMatrix();
}

void chair()
{
    glPushMatrix();
    glTranslatef(1,0,0.3);

    glPushMatrix();
    glTranslatef(6,0,3);
    glScalef(0.07,1.2,.07);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(6.6,0,3);
    glScalef(0.07,1.2,.07);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(6.6,0,3.6);
    glScalef(0.07,2.2,.07);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(6,0,3.6);
    glScalef(0.07,2.2,.07);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(5.9,1.2,3);
    glScalef(.8,.07,.8);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(5.8,2.2,3.6);
    glScalef(1,1,.07);
    cube(1.000, 1.000, 0.000);
    glPopMatrix();


    glPopMatrix();


    ///chair 2

    glPushMatrix();
    glTranslatef(-13.8,0,0.3);

    glPushMatrix();
    glTranslatef(6,0,3);
    glScalef(0.07,1.2,.07);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(6.6,0,3);
    glScalef(0.07,1.2,.07);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(6.6,0,3.6);
    glScalef(0.07,2.2,.07);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(6,0,3.6);
    glScalef(0.07,2.2,.07);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(5.9,1.2,3);
    glScalef(.8,.07,.8);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(5.8,2.2,3.6);
    glScalef(1,1,.07);
    cube(1.000, 1.000, 0.000);
    glPopMatrix();


    glPopMatrix();


}

void field(GLfloat l,GLfloat w)
{
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,1);
    glPushMatrix();
    glTranslatef(-(l/2),0,-(w/2));
    glScalef(l,0.1,w);
    cube(1,1,1);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);
}
void umbrella()
{


    set_mat_prop(0.596, 0.984, 0.596,true,60);


    if(wired)
    {
        glPolygonMode( GL_FRONT, GL_LINE ) ;
        glPolygonMode( GL_BACK, GL_LINE ) ;

    }
    else
    {
        glPolygonMode( GL_FRONT,GL_FILL ) ;
        glPolygonMode( GL_BACK, GL_FILL ) ;
    }





    glPushMatrix();
    glPushMatrix();
    glTranslatef(4,8,-5);
    glRotated(-90,0,0,1);
    glColor3f(2,2,2);
    bottleBezier();
    glPopMatrix();

    glPushMatrix();
    glTranslatef(3.5,0,-3.5);
    glScalef(.2,7,.1);
    cube(0.333, 0.420, 0.184);
    glPopMatrix();

    glPopMatrix();

    glPushMatrix();
    set_mat_prop(0.000, 1.000, 0.498,true,60);
    glTranslatef(0,17.5,-2);
    glRotated(-90,0,0,1);
    glScaled(5.5,5.5,5.5);

    bottleBezier();
    glPopMatrix();


    ///food court
    glPushMatrix();

    glPushMatrix();
    glTranslated(6,0,18);
    glPushMatrix();
    glTranslatef(4,8,-5);
    glRotated(-90,0,0,1);
    glColor3f(2,2,2);
    bottleBezier();
    glPopMatrix();

    glPushMatrix();
    glTranslatef(3.5,0,-3.5);
    glScalef(.2,7,.1);
    cube(0.333, 0.420, 0.184);
    glPopMatrix();

    glPopMatrix();
    ///table

    glPushMatrix();
    glTranslated(8.5,2,13.5);
    glScaled(2,.2,2);
    cube(0.604, 0.804, 0.196);

    glPopMatrix();

    ///chair
    glPushMatrix();
    glTranslatef(3.5,0,12);

    glPushMatrix();
    glTranslatef(6,0,3);
    glScalef(0.07,1.2,.07);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(6.6,0,3);
    glScalef(0.07,1.2,.07);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(6.6,0,3.6);
    glScalef(0.07,2.2,.07);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(6,0,3.6);
    glScalef(0.07,2.2,.07);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(5.9,1.2,3);
    glScalef(.8,.07,.8);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(5.8,2.2,3.6);
    glScalef(1,1,.07);
    cube(1.000, 1.000, 0.000);
    glPopMatrix();


    glPopMatrix();

    ///chair 2

    glPushMatrix();

    glTranslatef(16,0,16);
    glRotated(180,0,1,0);


    glPushMatrix();
    glTranslatef(6,0,3);
    glScalef(0.07,1.2,.07);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(6.6,0,3);
    glScalef(0.07,1.2,.07);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(6.6,0,3.6);
    glScalef(0.07,2.2,.07);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(6,0,3.6);
    glScalef(0.07,2.2,.07);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(5.9,1.2,3);
    glScalef(.8,.07,.8);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(5.8,2.2,3.6);
    glScalef(1,1,.07);
    cube(1.000, 1.000, 0.000);
    glPopMatrix();


    glPopMatrix();

    glPopMatrix();


     ///food court 2
    glPushMatrix();
    glTranslated(-16,0,0);
    glPushMatrix();
    set_mat_prop(0.804, 0.522, 0.247,true,60);
    glTranslated(6,0,18);
    glPushMatrix();
    glTranslatef(4,8,-5);
    glRotated(-90,0,0,1);
    glColor3f(2,2,2);
    bottleBezier();
    glPopMatrix();

    glPushMatrix();
    glTranslatef(3.5,0,-3.5);
    glScalef(.2,7,.1);
    cube(0.333, 0.420, 0.184);
    glPopMatrix();

    glPopMatrix();
    ///table

    glPushMatrix();
    glTranslated(8.5,2,13.5);
    glScaled(2,.2,2);
    cube(0.737, 0.561, 0.561);

    glPopMatrix();

    ///chair
    glPushMatrix();
    glTranslatef(3.5,0,12);

    glPushMatrix();
    glTranslatef(6,0,3);
    glScalef(0.07,1.2,.07);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(6.6,0,3);
    glScalef(0.07,1.2,.07);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(6.6,0,3.6);
    glScalef(0.07,2.2,.07);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(6,0,3.6);
    glScalef(0.07,2.2,.07);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(5.9,1.2,3);
    glScalef(.8,.07,.8);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(5.8,2.2,3.6);
    glScalef(1,1,.07);
    cube(1.000, 1.000, 0.000);
    glPopMatrix();


    glPopMatrix();

    ///chair 2

    glPushMatrix();

    glTranslatef(16,0,16);
    glRotated(180,0,1,0);


    glPushMatrix();
    glTranslatef(6,0,3);
    glScalef(0.07,1.2,.07);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(6.6,0,3);
    glScalef(0.07,1.2,.07);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(6.6,0,3.6);
    glScalef(0.07,2.2,.07);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(6,0,3.6);
    glScalef(0.07,2.2,.07);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(5.9,1.2,3);
    glScalef(.8,.07,.8);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(5.8,2.2,3.6);
    glScalef(1,1,.07);
    cube(1.000, 1.000, 0.000);
    glPopMatrix();


    glPopMatrix();

    glPopMatrix();


    ///food court 3

    glPushMatrix();
    glTranslated(-17,0,7.5);
    glPushMatrix();
    set_mat_prop(0.957, 0.643, 0.376,true,60);
    glTranslated(6,0,18);
    glPushMatrix();
    glTranslatef(4,8,-5);
    glRotated(-90,0,0,1);
    glColor3f(2,2,2);
    bottleBezier();
    glPopMatrix();

    glPushMatrix();
    glTranslatef(3.5,0,-3.5);
    glScalef(.2,7,.1);
    cube(0.333, 0.420, 0.184);
    glPopMatrix();

    glPopMatrix();
    ///table

    glPushMatrix();
    glTranslated(8.5,2,13.5);
    glScaled(2,.2,2);
    cube(0.855, 0.647, 0.125);

    glPopMatrix();

    ///chair
    glPushMatrix();
    glTranslatef(3.5,0,12);

    glPushMatrix();
    glTranslatef(6,0,3);
    glScalef(0.07,1.2,.07);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(6.6,0,3);
    glScalef(0.07,1.2,.07);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(6.6,0,3.6);
    glScalef(0.07,2.2,.07);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(6,0,3.6);
    glScalef(0.07,2.2,.07);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(5.9,1.2,3);
    glScalef(.8,.07,.8);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(5.8,2.2,3.6);
    glScalef(1,1,.07);
    cube(1.000, 1.000, 0.000);
    glPopMatrix();


    glPopMatrix();

    ///chair 2

    glPushMatrix();

    glTranslatef(16,0,16);
    glRotated(180,0,1,0);


    glPushMatrix();
    glTranslatef(6,0,3);
    glScalef(0.07,1.2,.07);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(6.6,0,3);
    glScalef(0.07,1.2,.07);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(6.6,0,3.6);
    glScalef(0.07,2.2,.07);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(6,0,3.6);
    glScalef(0.07,2.2,.07);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(5.9,1.2,3);
    glScalef(.8,.07,.8);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(5.8,2.2,3.6);
    glScalef(1,1,.07);
    cube(1.000, 1.000, 0.000);
    glPopMatrix();


    glPopMatrix();

    glPopMatrix();

    ///food court 4



    glPushMatrix();
    glTranslated(0,0,7.5);
    glPushMatrix();
    set_mat_prop(0.941, 1.000, 0.941,true,60);
    glTranslated(6,0,18);
    glPushMatrix();
    glTranslatef(4,8,-5);
    glRotated(-90,0,0,1);
    glColor3f(2,2,2);
    bottleBezier();
    glPopMatrix();

    glPushMatrix();
    glTranslatef(3.5,0,-3.5);
    glScalef(.2,7,.1);
    cube(0.333, 0.420, 0.184);
    glPopMatrix();

    glPopMatrix();
    ///table

    glPushMatrix();
    glTranslated(8.5,2,13.5);
    glScaled(2,.2,2);
    cube(	0.467, 0.533, 0.600);

    glPopMatrix();

    ///chair
    glPushMatrix();
    glTranslatef(3.5,0,12);

    glPushMatrix();
    glTranslatef(6,0,3);
    glScalef(0.07,1.2,.07);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(6.6,0,3);
    glScalef(0.07,1.2,.07);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(6.6,0,3.6);
    glScalef(0.07,2.2,.07);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(6,0,3.6);
    glScalef(0.07,2.2,.07);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(5.9,1.2,3);
    glScalef(.8,.07,.8);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(5.8,2.2,3.6);
    glScalef(1,1,.07);
    cube(1.000, 1.000, 0.000);
    glPopMatrix();


    glPopMatrix();

    ///chair 2

    glPushMatrix();

    glTranslatef(16,0,16);
    glRotated(180,0,1,0);


    glPushMatrix();
    glTranslatef(6,0,3);
    glScalef(0.07,1.2,.07);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(6.6,0,3);
    glScalef(0.07,1.2,.07);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(6.6,0,3.6);
    glScalef(0.07,2.2,.07);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(6,0,3.6);
    glScalef(0.07,2.2,.07);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(5.9,1.2,3);
    glScalef(.8,.07,.8);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(5.8,2.2,3.6);
    glScalef(1,1,.07);
    cube(1.000, 1.000, 0.000);
    glPopMatrix();


    glPopMatrix();

    glPopMatrix();







    if(shcpt)
    {

        showControlPoints();
    }



}
void road()
{


    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,2);
    glPushMatrix();
    glTranslatef(-12,0,-12);
    glScalef(4,0.1,17);
    cube(1,1,1);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,2);
    glPushMatrix();
    glTranslatef(-12,0,4);
    glScalef(24,0.1,3);
    cube(1,1,1);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);


    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,2);
    glPushMatrix();
    glTranslatef(8,0,-12);
    glScalef(4,0.1,17);
    cube(1,1,1);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,2);
    glPushMatrix();
    glTranslatef(-3,0,7);
    glScalef(8.5,0.1,17);
    cube(1,1,1);
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);



    glPushMatrix();

    glTranslated(1,0,7);
    glScaled(.4,.2,17);
    cube(1,1,1);
    glPopMatrix();
}



void gallery(GLfloat l,GLfloat w)
{
    glPushMatrix();
    glTranslatef(0,0,-1);
    glPushMatrix();
    glTranslatef(-(l/2)-0.125,0,-(w/2));
    glScalef(l,1,.8);
    cube(1,1,1);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-(l/2)-0.125,1,-(w/2)-.8);
    glScalef(l,1,.8);
    cube(0.980, 0.502, 0.447);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-(l/2)-0.125,2,-(w/2)-2*.8);
    glScalef(l,1,.8);
    cube(1,1,1);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-(l/2)-0.125,3,-(w/2)-3*.8);
    glScalef(l,1,.8);
    cube(0.980, 0.502, 0.447);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-(l/2)-0.125,4,-(w/2)-4*.8);
    glScalef(l,1,.8);
    cube(1,1,1);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-(l/2)-0.125,5,-(w/2)-5*.8);
    glScalef(l,1,.8);
    cube(0.980, 0.502, 0.447);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-(l/2)-0.125,6,-(w/2)-6*.8);
    glScalef(l,1,.8);
    cube(1,1,1);
    glPopMatrix();

    glPopMatrix();

}
void tree()
{
    set_mat_prop(0,.4,0,true,60);



    glPushMatrix();
    glTranslated(-6.7,3,3.2);
    glPushMatrix();

    ///glTranslated(-6.5,3,3.5);
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    GLUquadricObj *quadratic;
    quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.125,0.125,10,20,8);
    glPopMatrix();

    glPushMatrix();
    glTranslated(0,.5,0);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();



    glPushMatrix();
    glTranslated(0,-.5,0);
    glRotated(-45,0,0,1);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();


    glPushMatrix();
    glTranslated(0,-.5,0);
    glRotated(45,0,0,1);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();
    ///
    glPopMatrix();


    ///
    glPopMatrix();



    ///2nd tree

    glPushMatrix();
    glTranslated(-6.7,3,-1.5);
    glPushMatrix();

    ///glTranslated(-6.5,3,3.5);
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    /// quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.125,0.125,10,20,8);
    glPopMatrix();

    glPushMatrix();
    glTranslated(0,.5,0);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();



    glPushMatrix();
    glTranslated(0,-.5,0);
    glRotated(-45,0,0,1);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();


    glPushMatrix();
    glTranslated(0,-.5,0);
    glRotated(45,0,0,1);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();
    ///
    glPopMatrix();


    ///
    glPopMatrix();


///3rd tree

    glPushMatrix();
    glTranslated(8,3,-1.5);
    glPushMatrix();

    ///glTranslated(-6.5,3,3.5);
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    /// quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.125,0.125,10,20,8);
    glPopMatrix();

    glPushMatrix();
    glTranslated(0,.5,0);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();



    glPushMatrix();
    glTranslated(0,-.5,0);
    glRotated(-45,0,0,1);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();


    glPushMatrix();
    glTranslated(0,-.5,0);
    glRotated(45,0,0,1);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();
    ///
    glPopMatrix();


    ///
    glPopMatrix();


///4th tree

    glPushMatrix();
    glTranslated(8,3,1);
    glPushMatrix();

    ///glTranslated(-6.5,3,3.5);
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    /// quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.125,0.125,10,20,8);
    glPopMatrix();

    glPushMatrix();
    glTranslated(0,.5,0);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();



    glPushMatrix();
    glTranslated(0,-.5,0);
    glRotated(-45,0,0,1);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();


    glPushMatrix();
    glTranslated(0,-.5,0);
    glRotated(45,0,0,1);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();
    ///
    glPopMatrix();


    ///
    glPopMatrix();

    ///5th tree

    glPushMatrix();
    glTranslated(8,3,3);
    glPushMatrix();

    ///glTranslated(-6.5,3,3.5);
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    /// quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.125,0.125,10,20,8);
    glPopMatrix();

    glPushMatrix();
    glTranslated(0,.5,0);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();



    glPushMatrix();
    glTranslated(0,-.5,0);
    glRotated(-45,0,0,1);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();


    glPushMatrix();
    glTranslated(0,-.5,0);
    glRotated(45,0,0,1);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();
    ///
    glPopMatrix();


    ///
    glPopMatrix();




    ///road side tree 1

    glPushMatrix();
    glTranslated(-2.2,3,10);
    glPushMatrix();

    ///glTranslated(-6.5,3,3.5);
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    /// quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.125,0.125,10,20,8);
    glPopMatrix();

    glPushMatrix();
    glTranslated(0,.5,0);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();



    glPushMatrix();
    glTranslated(0,-.5,0);
    glRotated(-45,0,0,1);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();


    glPushMatrix();
    glTranslated(0,-.5,0);
    glRotated(45,0,0,1);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();
    ///
    glPopMatrix();


    ///
    glPopMatrix();



    ///road side tree 2

    glPushMatrix();
    glTranslated(-2.2,3,14);
    glPushMatrix();

    ///glTranslated(-6.5,3,3.5);
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    /// quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.125,0.125,10,20,8);
    glPopMatrix();

    glPushMatrix();
    glTranslated(0,.5,0);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();



    glPushMatrix();
    glTranslated(0,-.5,0);
    glRotated(-45,0,0,1);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();


    glPushMatrix();
    glTranslated(0,-.5,0);
    glRotated(45,0,0,1);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();
    ///
    glPopMatrix();


    ///
    glPopMatrix();

    ///road side tree 3

    glPushMatrix();
    glTranslated(-2.2,3,18);
    glPushMatrix();

    ///glTranslated(-6.5,3,3.5);
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    /// quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.125,0.125,10,20,8);
    glPopMatrix();

    glPushMatrix();
    glTranslated(0,.5,0);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();



    glPushMatrix();
    glTranslated(0,-.5,0);
    glRotated(-45,0,0,1);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();


    glPushMatrix();
    glTranslated(0,-.5,0);
    glRotated(45,0,0,1);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();
    ///
    glPopMatrix();


    ///
    glPopMatrix();

    ///road side tree 4

    glPushMatrix();
    glTranslated(-2.2,3,22);
    glPushMatrix();

    ///glTranslated(-6.5,3,3.5);
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    /// quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.125,0.125,10,20,8);
    glPopMatrix();

    glPushMatrix();
    glTranslated(0,.5,0);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();



    glPushMatrix();
    glTranslated(0,-.5,0);
    glRotated(-45,0,0,1);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();


    glPushMatrix();
    glTranslated(0,-.5,0);
    glRotated(45,0,0,1);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();
    ///
    glPopMatrix();


    ///
    glPopMatrix();


    ///road side tree 5

    glPushMatrix();
    glTranslated(-2.2,3,25);
    glPushMatrix();

    ///glTranslated(-6.5,3,3.5);
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    /// quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.125,0.125,10,20,8);
    glPopMatrix();

    glPushMatrix();
    glTranslated(0,.5,0);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();



    glPushMatrix();
    glTranslated(0,-.5,0);
    glRotated(-45,0,0,1);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();


    glPushMatrix();
    glTranslated(0,-.5,0);
    glRotated(45,0,0,1);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();
    ///
    glPopMatrix();


    ///
    glPopMatrix();


    ///right side tree

    ///glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();
    ///
    glPopMatrix();


    ///
    glPopMatrix();

    ///road side tree 6

    glPushMatrix();
    set_mat_prop(0.000, 1.000, 0.000,true,60);
    glTranslated(-2.2,3,25);
    glPushMatrix();

    ///glTranslated(-6.5,3,3.5);
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    /// quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.125,0.125,10,20,8);
    glPopMatrix();

    glPushMatrix();
    glTranslated(0,.5,0);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();



    glPushMatrix();
    glTranslated(0,-.5,0);
    glRotated(-45,0,0,1);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();


    glPushMatrix();
    glTranslated(0,-.5,0);
    glRotated(45,0,0,1);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();
    ///
    glPopMatrix();


    ///
    glPopMatrix();

    ///glPopMatrix();

    ///right side tree  1

    glPushMatrix();
    set_mat_prop(0.000, 1.000, 0.000,true,60);
    glTranslated(4.2,3,22);
    glPushMatrix();

    ///glTranslated(-6.5,3,3.5);
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    /// quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.125,0.125,10,20,8);
    glPopMatrix();

    glPushMatrix();
    glTranslated(0,.5,0);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();



    glPushMatrix();
    glTranslated(0,-.5,0);
    glRotated(-45,0,0,1);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();


    glPushMatrix();
    glTranslated(0,-.5,0);
    glRotated(45,0,0,1);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();
    ///
    glPopMatrix();


    ///
    glPopMatrix();



    ///right side tree  2

    glPushMatrix();
    set_mat_prop(0.000, 1.000, 0.000,true,60);
    glTranslated(4.2,3,14);
    glPushMatrix();

    ///glTranslated(-6.5,3,3.5);
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    /// quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.125,0.125,10,20,8);
    glPopMatrix();

    glPushMatrix();
    glTranslated(0,.5,0);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();



    glPushMatrix();
    glTranslated(0,-.5,0);
    glRotated(-45,0,0,1);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();


    glPushMatrix();
    glTranslated(0,-.5,0);
    glRotated(45,0,0,1);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();
    ///
    glPopMatrix();


    ///
    glPopMatrix();



    ///right side tree  3

    glPushMatrix();
    glTranslated(4.2,3,18);
    glPushMatrix();

    ///glTranslated(-6.5,3,3.5);
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    /// quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.125,0.125,10,20,8);
    glPopMatrix();

    glPushMatrix();
    glTranslated(0,.5,0);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();



    glPushMatrix();
    glTranslated(0,-.5,0);
    glRotated(-45,0,0,1);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();


    glPushMatrix();
    glTranslated(0,-.5,0);
    glRotated(45,0,0,1);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();
    ///
    glPopMatrix();


    ///
    glPopMatrix();


    ///right side tree  1

    glPushMatrix();
    set_mat_prop(0.000, 1.000, 0.000,true,60);
    glTranslated(4.2,3,10);
    glPushMatrix();

    ///glTranslated(-6.5,3,3.5);
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    /// quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.125,0.125,10,20,8);
    glPopMatrix();

    glPushMatrix();
    glTranslated(0,.5,0);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();



    glPushMatrix();
    glTranslated(0,-.5,0);
    glRotated(-45,0,0,1);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();


    glPushMatrix();
    glTranslated(0,-.5,0);
    glRotated(45,0,0,1);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();
    ///
    glPopMatrix();


    ///
    glPopMatrix();




    ///right side tree  5

    glPushMatrix();
    set_mat_prop(0.000, 1.000, 0.000,true,60);
    glTranslated(4.2,3,25);
    glPushMatrix();

    ///glTranslated(-6.5,3,3.5);
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    /// quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.125,0.125,10,20,8);
    glPopMatrix();

    glPushMatrix();
    glTranslated(0,.5,0);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();



    glPushMatrix();
    glTranslated(0,-.5,0);
    glRotated(-45,0,0,1);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();


    glPushMatrix();
    glTranslated(0,-.5,0);
    glRotated(45,0,0,1);
    ///daty
    glPushMatrix();
    glScaled(1,0.2,1);
    glRotated(90,1,0,0);

    ///GLUquadricObj *quadratic;
    ///quadratic=gluNewQuadric();
    gluCylinder(quadratic,0.0125,0.0125,3,20,8);


    glPopMatrix();
    ///leave
    glPushMatrix();
    glTranslated(0,.5,0);
    glutSolidSphere(0.5,20,20);

    glPopMatrix();
    ///
    glPopMatrix();
    ///
    glPopMatrix();


    ///
    glPopMatrix();











}
void car()
{
    glPushMatrix();
    glTranslated(-.25,0,c);
    glPushMatrix();
    glTranslated(0,0,18);
    glScaled(3,1.5,4);
    cube(0.627, 0.322, 0.176);
    glPopMatrix();

    glPushMatrix();
    glTranslated(0,1.5,18);
    glScaled(3,2,1.5);
    cube(1,1,0);
    glPopMatrix();

    glPushMatrix();
    glTranslated(0,1.5,19.5);
    glScaled(.2,1.4,2.5);
    cube(0.000, 0.502, 0.502);
    glPopMatrix();

    glPushMatrix();
    glTranslated(2.8,1.5,19.5);
    glScaled(.2,1.4,2.5);
    cube(0.000, 0.502, 0.502);
    glPopMatrix();

    ///car 2






     ///bottle

    glPushMatrix();
    set_mat_prop(1,0,0,true,60);
    glTranslatef(1.5,1.5,21);
    glRotated(90,0,0,1);
    glScaled(.9,.8,1);

    bottleBezier1();
    glPopMatrix();


    glPopMatrix();





}

void banner(GLfloat l,GLfloat w)
{
    glPushMatrix();


    ///left wall

    /*glPushMatrix();
    glTranslatef(-(l/2)-0.125,0,-(w/2));
    glScalef(0.125,1.5,w);
    cube(0.486, 0.988, 0.000);
    glPopMatrix();*/

    ///right wall

    glPushMatrix();
    glTranslatef(8,0,-(w/2));
    glScalef(0.125,1.5,w);
    cube(0.486, 0.988, 0.000);
    glPopMatrix();
    glPopMatrix();
}

void lamp(GLfloat l,GLfloat w)
{
    glPushMatrix();


    ///right top


    glPushMatrix();
    glTranslatef(8-.125,0,-(w/2));
    glScalef(.2,7,.2);
    cube(0.871, 0.722, 0.529);
    glPopMatrix();

    ///light board1
    glPushMatrix();
    glTranslatef(6.5,7,-(w/2));
    glScalef(3,2.5,.125);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();



    ///left top

    glPushMatrix();
    glTranslatef(-(l/2)-0.125,0,-(w/2));
    glScalef(.2,7,.2);
    cube(0.871, 0.722, 0.529);
    glPopMatrix();

    ///light board2
    glPushMatrix();
    glTranslatef(-9.5,7,-(w/2));
    glScalef(3,2.5,.125);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();


    ///left bottom

    glPushMatrix();
    glTranslatef(-(l/2)-0.125,0,(w/2));
    glScalef(.2,6,.2);
    cube(0.871, 0.722, 0.529);
    glPopMatrix();

    ///light board3
    glPushMatrix();
    glTranslatef(-9.5,6,(w/2));
    glScalef(3,2.5,.125);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();

    ///right bottom

    glPushMatrix();
    glTranslatef(8-.125,0,(w/2));
    glScalef(.2,6,.2);
    cube(0.871, 0.722, 0.529);
    glPopMatrix();

    ///light board4
    glPushMatrix();
    glTranslatef(6.5,6,(w/2));
    glScalef(3,2.5,.125);
    cube(0.502, 0.502, 0.000);
    glPopMatrix();



    glPopMatrix();

}
void player()
{
    /* glPushMatrix();
     glScaled(.8,.5,1);


     glPushMatrix();
     glTranslated(0,2.5,2);
     glScaled(1.5,2,.1);
     cube(0.545, 0.271, 0.075);
     glPopMatrix();

     glPushMatrix();
     glTranslated(0.2,0,2);
     glScaled(.2,2.5,.1);
     cube(0.545, 0.271, 0.075);
     glPopMatrix();

     glPushMatrix();
     glTranslated(1.1,0,2);
     glScaled(.2,2.5,.1);
     cube(0.545, 0.271, 0.075);
     glPopMatrix();

     glPushMatrix();
     set_mat_prop(1,1,1,true,60);
     glTranslated(.8,5.1,2);
     glutSolidSphere(0.6,20,20);
     glPopMatrix();

     glPushMatrix();
     glTranslated(1.5,4,2);
     glScaled(1,.3,.1);
     cube(0.545, 0.271, 0.075);
     glPopMatrix();

     glPushMatrix();
     glTranslated(-1,4,2);
     glScaled(1,.3,.1);
     cube(0.545, 0.271, 0.075);
     glPopMatrix();


     glPushMatrix();
     glTranslated(1,5.5,2.57);
     glScaled(.1,.1,.1);
     cube(0,0,0);
     glPopMatrix();

     glPushMatrix();
     glTranslated(.5,5.5,2.57);
     glScaled(.1,.1,.1);
     cube(0,0,0);
     glPopMatrix();


     glPopMatrix();*/


    ///2nd man

    glPushMatrix();
    glTranslated(6,0,4);

    glRotated(-135,0,1,0);
    glScaled(.8,.5,1);


    glPushMatrix();
    glTranslated(0,2.5,2);
    glScaled(1.5,2,.1);
    cube(0.545, 0.271, 0.075);
    glPopMatrix();

    glPushMatrix();
    glTranslated(0.2,0,2);
    glScaled(.2,2.5,.1);
    cube(0.545, 0.271, 0.075);
    glPopMatrix();

    glPushMatrix();
    glTranslated(1.1,0,2);
    glScaled(.2,2.5,.1);
    cube(0.545, 0.271, 0.075);
    glPopMatrix();

    glPushMatrix();
    set_mat_prop(1,1,1,true,60);
    glTranslated(.8,5.1,2);
    glutSolidSphere(0.6,20,20);
    glPopMatrix();

    glPushMatrix();
    glTranslated(1.5,4,2);
    glScaled(1,.3,.1);
    cube(0.545, 0.271, 0.075);
    glPopMatrix();

    glPushMatrix();
    glTranslated(-1,4,2);
    glScaled(1,.3,.1);
    cube(0.545, 0.271, 0.075);
    glPopMatrix();


    glPushMatrix();
    glTranslated(1,5.5,2.3);
    glScaled(.1,.1,.1);
    cube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslated(.5,5.5,2.3);
    glScaled(.1,.1,.1);
    cube(0,0,0);
    glPopMatrix();


    glPopMatrix();


    ///3rd man

    glPushMatrix();
    glTranslated(4,0,-3);
    glRotated(-20,0,1,0);
    glScaled(.8,.5,1);



    glPushMatrix();
    glTranslated(0,2.5,2);
    glScaled(1.5,2,.1);
    cube(0.545, 0.271, 0.075);
    glPopMatrix();

    glPushMatrix();
    glTranslated(0.2,0,2);
    glScaled(.2,2.5,.1);
    cube(0.545, 0.271, 0.075);
    glPopMatrix();

    glPushMatrix();
    glTranslated(1.1,0,2);
    glScaled(.2,2.5,.1);
    cube(0.545, 0.271, 0.075);
    glPopMatrix();

    glPushMatrix();
    set_mat_prop(1,1,1,true,60);
    glTranslated(.8,5.1,2);
    glutSolidSphere(0.6,20,20);
    glPopMatrix();

    glPushMatrix();
    glTranslated(1.5,4,2);
    glScaled(1,.3,.1);
    cube(0.545, 0.271, 0.075);
    glPopMatrix();

    glPushMatrix();
    glTranslated(-1,4,2);
    glScaled(1,.3,.1);
    cube(0.545, 0.271, 0.075);
    glPopMatrix();


    glPushMatrix();
    glTranslated(1,5.5,2.3);
    glScaled(.1,.1,.1);
    cube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslated(.5,5.5,2.3);
    glScaled(.1,.1,.1);
    cube(0,0,0);
    glPopMatrix();


    glPopMatrix();


    ///4th man



    ///4th man

    /* glPushMatrix();
     glScaled(.8,.5,1);
     glTranslated(0,0,-3);

     glPushMatrix();
     glTranslated(0,2.5,2);
     glScaled(1.5,2,.1);
     cube(0.545, 0.271, 0.075);
     glPopMatrix();

     glPushMatrix();
     glTranslated(0.2,0,2);
     glScaled(.2,2.5,.1);
     cube(0.545, 0.271, 0.075);
     glPopMatrix();

     glPushMatrix();
     glTranslated(1.1,0,2);
     glScaled(.2,2.5,.1);
     cube(0.545, 0.271, 0.075);
     glPopMatrix();

     glPushMatrix();
     set_mat_prop(1,1,1,true,60);
     glTranslated(.8,5.1,2);
     glutSolidSphere(0.6,20,20);
     glPopMatrix();

     glPushMatrix();
     glTranslated(1.5,4,2);
     glScaled(1,.3,.1);
     cube(0.545, 0.271, 0.075);
     glPopMatrix();

     glPushMatrix();
     glTranslated(-1,4,2);
     glScaled(1,.3,.1);
     cube(0.545, 0.271, 0.075);
     glPopMatrix();


     glPushMatrix();
     glTranslated(1,5.5,2.57);
     glScaled(.1,.1,.1);
     cube(0,0,0);
     glPopMatrix();

     glPushMatrix();
     glTranslated(.5,5.5,2.57);
     glScaled(.1,.1,.1);
     cube(0,0,0);
     glPopMatrix();


     glPopMatrix();*/


    ///5th man

    glPushMatrix();

    glTranslated(-5,0,-3);
    glRotated(20,0,1,0);

    glScaled(.8,.5,1);

    glPushMatrix();
    glTranslated(0,2.5,2);
    glScaled(1.5,2,.1);
    cube(0.855, 0.647, 0.125);
    glPopMatrix();

    glPushMatrix();
    glTranslated(0.2,0,2);
    glScaled(.2,2.5,.1);
    cube(0.855, 0.647, 0.125);
    glPopMatrix();

    glPushMatrix();
    glTranslated(1.1,0,2);
    glScaled(.2,2.5,.1);
    cube(0.855, 0.647, 0.125);
    glPopMatrix();

    glPushMatrix();
    set_mat_prop(1,1,1,true,60);
    glTranslated(.8,5.1,2);
    glutSolidSphere(0.6,20,20);
    glPopMatrix();

    glPushMatrix();
    glTranslated(1.5,4,2);
    glScaled(1,.3,.1);
    cube(0.855, 0.647, 0.125);
    glPopMatrix();

    glPushMatrix();
    glTranslated(-1,4,2);
    glScaled(1,.3,.1);
    cube(0.855, 0.647, 0.125);
    glPopMatrix();


    glPushMatrix();
    glTranslated(1,5.5,2.3);
    glScaled(.1,.1,.1);
    cube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslated(.5,5.5,2.3);
    glScaled(.1,.1,.1);
    cube(0,0,0);
    glPopMatrix();


    glPopMatrix();


    ///6th man

    /* glPushMatrix();
    glScaled(.8,.5,1);
    glTranslated(-7.5,0,-3);

    glPushMatrix();
    glTranslated(0,2.5,2);
    glScaled(1.5,2,.1);
    cube(0.855, 0.647, 0.125);
    glPopMatrix();

    glPushMatrix();
    glTranslated(0.2,0,2);
    glScaled(.2,2.5,.1);
    cube(0.855, 0.647, 0.125);
    glPopMatrix();

    glPushMatrix();
    glTranslated(1.1,0,2);
    glScaled(.2,2.5,.1);
    cube(0.855, 0.647, 0.125);
    glPopMatrix();

    glPushMatrix();
    set_mat_prop(1,1,1,true,60);
    glTranslated(.8,5.1,2);
    glutSolidSphere(0.6,20,20);
    glPopMatrix();

    glPushMatrix();
    glTranslated(1.5,4,2);
    glScaled(1,.3,.1);
    cube(0.855, 0.647, 0.125);
    glPopMatrix();

    glPushMatrix();
    glTranslated(-1,4,2);
    glScaled(1,.3,.1);
    cube(0.855, 0.647, 0.125);
    glPopMatrix();


    glPushMatrix();
    glTranslated(1,5.5,2.57);
    glScaled(.1,.1,.1);
    cube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslated(.5,5.5,2.57);
    glScaled(.1,.1,.1);
    cube(0,0,0);
    glPopMatrix();


    glPopMatrix();*/


    ///7th man


    /*glPushMatrix();
    glScaled(.8,.5,1);
    glTranslated(-5,0,0);

    glPushMatrix();
    glTranslated(0,2.5,2);
    glScaled(1.5,2,.1);
    cube(0.855, 0.647, 0.125);
    glPopMatrix();

    glPushMatrix();
    glTranslated(0.2,0,2);
    glScaled(.2,2.5,.1);
    cube(0.855, 0.647, 0.125);
    glPopMatrix();

    glPushMatrix();
    glTranslated(1.1,0,2);
    glScaled(.2,2.5,.1);
    cube(0.855, 0.647, 0.125);
    glPopMatrix();

    glPushMatrix();
    set_mat_prop(1,1,1,true,60);
    glTranslated(.8,5.1,2);
    glutSolidSphere(0.6,20,20);
    glPopMatrix();

    glPushMatrix();
    glTranslated(1.5,4,2);
    glScaled(1,.3,.1);
    cube(0.855, 0.647, 0.125);
    glPopMatrix();

    glPushMatrix();
    glTranslated(-1,4,2);
    glScaled(1,.3,.1);
    cube(0.855, 0.647, 0.125);
    glPopMatrix();


    glPushMatrix();
    glTranslated(1,5.5,2.3);
    glScaled(.1,.1,.1);
    cube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslated(.5,5.5,2.3);
    glScaled(.1,.1,.1);
    cube(0,0,0);
    glPopMatrix();


    glPopMatrix();

    ///8th man

    glPushMatrix();
    glScaled(.8,.5,1);
    glTranslated(-3.7,0,0);

    glPushMatrix();
    glTranslated(0,2.5,2);
    glScaled(1.5,2,.1);
    cube(0.855, 0.647, 0.125);
    glPopMatrix();

    glPushMatrix();
    glTranslated(0.2,0,2);
    glScaled(.2,2.5,.1);
    cube(0.855, 0.647, 0.125);
    glPopMatrix();

    glPushMatrix();
    glTranslated(1.1,0,2);
    glScaled(.2,2.5,.1);
    cube(0.855, 0.647, 0.125);
    glPopMatrix();

    glPushMatrix();
    set_mat_prop(1,1,1,true,60);
    glTranslated(.8,5.1,2);
    glutSolidSphere(0.6,20,20);
    glPopMatrix();

    glPushMatrix();
    glTranslated(1.5,4,2);
    glScaled(1,.3,.1);
    cube(0.855, 0.647, 0.125);
    glPopMatrix();

    glPushMatrix();
    glTranslated(-1,4,2);
    glScaled(1,.3,.1);
    cube(0.855, 0.647, 0.125);
    glPopMatrix();


    glPushMatrix();
    glTranslated(1,5.5,2.57);
    glScaled(.1,.1,.1);
    cube(0,0,0);
    glPopMatrix();

    glPushMatrix();
    glTranslated(.5,5.5,2.57);
    glScaled(.1,.1,.1);
    cube(0,0,0);
    glPopMatrix();


    glPopMatrix(); */
















}
void goalbar()
{
/// right bar
    glPushMatrix();

    glPushMatrix();
    glTranslatef(7,0,-1);
    glScalef(0.1,3,.1);
    cube(1,1,1);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(7,0,1.5);
    glScalef(0.1,3,.1);
    cube(1,1,1);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(7,3,-.9);
    glScalef(0.1,.1,2.5);
    cube(1,1,1);
    glPopMatrix();

    glPopMatrix();

    /// left bar
    glPushMatrix();
    glTranslatef(-14,0,0);

    glPushMatrix();
    glTranslatef(7,0,-1);
    glScalef(0.1,3,.1);
    cube(1,1,1);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(7,0,1.5);
    glScalef(0.1,3,.1);
    cube(1,1,1);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(7,3,-.9);
    glScalef(0.1,.1,2.5);
    cube(1,1,1);
    glPopMatrix();

    glPopMatrix();



}

static void resize(int width, int height)
{
    float rat_new = 1.0*width/height;
    float width_new = height*rat;
    float height_new = width/rat;
    float x_t = 0;
    float y_t = 0;

    if(rat<rat_new)
    {
        x_t = (width-width_new)/2;
        width=width_new;
    }
    else
    {
        y_t = (height-height_new)/2;
        height=height_new;
    }

    glViewport(x_t, y_t, width, height);
}

static void display(void)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glFrustum(-4, 4, -3, 3, 3.0, 200.0);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(eyeX,eyeY,eyeZ, lookX,lookY,lookZ, 0,1,0);



    field(16,8);
    banner(16,8);
    lamp(16,8);
    gallery(16,8);
    goalbar();
    chair();
    umbrella();
    build();
    tree();
    road();
    player();
    ball();
    wall();
    car();
    showControlPoints();






    glLightfv(GL_LIGHT2, GL_AMBIENT,  l_amb2);
    glLightfv(GL_LIGHT2, GL_DIFFUSE,  l_dif2);
    glLightfv(GL_LIGHT2, GL_SPECULAR, l_spec2);
    glLightfv(GL_LIGHT2, GL_POSITION, l_pos2);
    //glLightfv(GL_LIGHT1, GL_POSITION, light_positions);
    if(spot_light==true)
    {
        glEnable(GL_LIGHT2);
        GLfloat l_spt[] = {0.0,0,-1,0.0};
        GLfloat spt_ct[] = {spot_cutoff};
        glLightfv(GL_LIGHT2, GL_SPOT_DIRECTION, l_spt);
        glLightfv(GL_LIGHT2, GL_SPOT_CUTOFF, spt_ct);
    }

    if(l_zero==true)
    {
        glEnable(GL_LIGHT0);

    }

    glFlush();
    glutSwapBuffers();
}

static void key(unsigned char key, int x, int y)
{
    float x_, z_, r, theta, dx, dz, dx_norm, dz_norm, r_=1, turn_angle_step=5, height_diff_one_less, height_diff_thresh_dist;

    x_=lookX-eyeX;
    z_=lookZ-eyeZ;
    r=sqrt(x_*x_+z_*z_);

    if(x_==0)
        theta = 90;
    else
        theta=atan(z_/x_) * 180 / PI;

    if((z_>0 && theta<0) || (z_<0 && theta>0))
        theta += 180;

    dx = r_*cos(theta * PI / 180);
    dz = r_*sin(theta * PI / 180);

    dx_norm = r_*cos((theta-90) * PI / 180);
    dz_norm = r_*sin((theta-90) * PI / 180);

    switch (key)
    {
    case 27 :
    case 'q':
        exit(0);
        break;

    case 'm':

        flagg=true;

        break;

    case 'M':
        flagg2=true;
        break;

    case 'n':
        flagg3=true;
        break;
    case 'N':
        flagg4=true;
        break;
    case 'v':
        flagg5=true;
        break;
    case 'V':
        flagg6=true;
        break;
    case 'x':
        flagg7=true;
        break;
    case 'X':
        flagg8=true;
        break;

    case 'u':
        flagg9=true;
        break;
    case 'U':
        flagg10=true;
        break;

    case 'o':
        spot_light=true;
        break;

    case 'O':
        spot_light=false;
        glDisable(GL_LIGHT2);
        break;

    // moving the look at point at a circular path or up-down
    case 'a':
        theta-=turn_angle_step;
        theta = theta * PI / 180;

        lookX=r*cos(theta)+eyeX;
        lookZ=r*sin(theta)+eyeZ;
        break;
    case 'w':
        lookY++;
        break;
    case 's':
        lookY--;
        break;
    case 'd':
        theta+=turn_angle_step;
        theta = theta * PI / 180;

        lookX=r*cos(theta)+eyeX;
        lookZ=r*sin(theta)+eyeZ;
        break;

    // Moving the camera front-back-left-right
    case 'j':
        eyeX += dx_norm;
        eyeZ += dz_norm;

        lookX += dx_norm;
        lookZ += dz_norm;
        break;
    case 'i':
        eyeX += dx;
        eyeZ += dz;

        lookX += dx;
        lookZ += dz;

        break;
    case 'k':
        eyeX -= dx;
        eyeZ -= dz;

        lookX -= dx;
        lookZ -= dz;
        break;
    case 'l':
        eyeX -= dx_norm;
        eyeZ -= dz_norm;

        lookX -= dx_norm;
        lookZ -= dz_norm;
        break;

    // Moving the camera up-down
    case '+':
        eyeY++;
        lookY++;
        break;
    case '-':
        eyeY--;
        lookY--;
        break;

    // Toggling window open-close
    case 't':
        wndw_rot = (((wndw_rot/89)+1) % 2) * 89;
        break;

    // rotating the whole scene
    case ',':
        rot--;
        break;
    case '.':
        rot++;
        break;

    // rotating the camera around the look at point
    case 'g':
        theta += 180;
        theta += turn_angle_step;
        theta = theta * PI / 180;

        eyeX=r*cos(theta)+lookX;
        eyeZ=r*sin(theta)+lookZ;
        break;
    case 'h':
        theta += 180;
        theta -= turn_angle_step;
        theta = theta * PI / 180;

        eyeX=r*cos(theta)+lookX;
        eyeZ=r*sin(theta)+lookZ;
        break;

    // look far or near
    case 'r':
        lookX += dx;
        lookZ += dz;
        break;
    case 'f':
        lookX -= dx;
        lookZ -= dz;
        break;

    // on-off switches to whichever it is attached


    // Ambient light on-off
    case '0':



       glEnable(GL_LIGHT0);

        glLightfv(GL_LIGHT0, GL_AMBIENT, l_amb);

        glLightfv(GL_LIGHT0, GL_DIFFUSE, l_dif);

        glLightfv(GL_LIGHT0, GL_SPECULAR, l_spec);

        glLightfv(GL_LIGHT0, GL_POSITION, l_pos);

        break;
    case '1':
        glDisable(GL_LIGHT0);
        glLightfv(GL_LIGHT0, GL_AMBIENT, l_no);





        break;

    case '2':
        glEnable(GL_LIGHT0);
        glLightfv(GL_LIGHT0,GL_AMBIENT,l_amb);
        glLightfv(GL_LIGHT0,GL_DIFFUSE,l_no);
        glLightfv(GL_LIGHT0,GL_SPECULAR,l_no);


        break;

    case '3':
        glEnable(GL_LIGHT0);
        glLightfv(GL_LIGHT0,GL_AMBIENT,l_no);
        glLightfv(GL_LIGHT0,GL_DIFFUSE,l_dif);
        glLightfv(GL_LIGHT0,GL_SPECULAR,l_no);


        break;

    case '4':
        glEnable(GL_LIGHT0);
        glLightfv(GL_LIGHT0,GL_AMBIENT,l_no);
        glLightfv(GL_LIGHT0,GL_DIFFUSE,l_no);
        glLightfv(GL_LIGHT0,GL_SPECULAR,l_spec);


        break;
    case '5':
        glEnable(GL_LIGHT1);

        glLightfv(GL_LIGHT1, GL_AMBIENT, l_amb1);

        glLightfv(GL_LIGHT1, GL_DIFFUSE, l_dif1);

        glLightfv(GL_LIGHT1, GL_SPECULAR, l_spec1);

        glLightfv(GL_LIGHT1, GL_POSITION, l_pos1);
        break;

    case '6':
        glDisable(GL_LIGHT1);
        break;

    case '7':
        glEnable(GL_LIGHT1);
        glLightfv(GL_LIGHT1,GL_AMBIENT,l_amb1);
        glLightfv(GL_LIGHT1,GL_DIFFUSE,l_no1);
        glLightfv(GL_LIGHT1,GL_SPECULAR,l_no1);
        break;
    case '8':
        glEnable(GL_LIGHT1);
        glLightfv(GL_LIGHT1,GL_AMBIENT,l_no1);
        glLightfv(GL_LIGHT1,GL_DIFFUSE,l_dif1);
        glLightfv(GL_LIGHT1,GL_SPECULAR,l_no1);
        break;

    case '9':
        glEnable(GL_LIGHT1);
        glLightfv(GL_LIGHT1,GL_AMBIENT,l_no1);
        glLightfv(GL_LIGHT1,GL_DIFFUSE,l_no1);
        glLightfv(GL_LIGHT1,GL_SPECULAR,l_spec1);
        break;



    case 'b':

        wired=!wired;
        break;


    }

    glutPostRedisplay();
}




int main(int argc, char *argv[])
{
    glutInit(&argc, argv);
    glutInitWindowSize(width,height);
    glutInitWindowPosition(30,30);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);

    glutCreateWindow("GLUT Shapes");
    glutIdleFunc(animate);

    glutDisplayFunc(display);
    glutKeyboardFunc(key);
    ///glutIdleFunc(idle);
    glutReshapeFunc(resize);

    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_SMOOTH);

    glEnable(GL_NORMALIZE);
    glEnable(GL_BLEND);
    glEnable(GL_LIGHTING);
    glutMouseFunc(processMouse);








    LoadTexture("C:\\Users\\USER\\Downloads\\yasin.bmp");
    LoadTexture("C:\\Users\\USER\\Downloads\\road.bmp");



    printf("Use 'w' to look up,\n 's' to look down,\n 'd' to look right,\n and 'a' to look left.\n\n");
    printf("Use 'i' to move camera front (i.e. zoom in),\n 'k' to move camera back (i.e. zoom out),\n 'l' to move camera right,\n and 'j' to move camera left.\n\n");
    printf("Use '+' to move camera up\n and '-' to move camera down.\n\n");
    printf("Use 't' to toggle window open or close.\n\n");
    printf("Use 'r' to look far,\n and 'f' to look near.\n\n");
    printf("Use 'g' to rotate left,\n and 'h' to rotate right taking the look at point as center.\n\n");
    printf("Use 'o' to toggle the target cross.\n\n");
    printf("Use 'F' to toggle the fan on-off.\n\n");
    printf("Use '1' and '2' to toggle tube lights switches,\n and '3' to toggle the table lamp switch.\n\n");
    printf("Use 'z' to toggle ambient,\n 'x' to toggle diffusion,\n and 'c' to toggle specular light property for all lights.\n\n\n\n");
    cout<<"Use 'm' for opening the door"<<endl;
    cout<<"Use 'M' for closing  the door"<<endl;
    cout<<"Use 'o' for  spot light"<<endl;
    cout<<"Use 'O' to stop the   spot light"<<endl;
    cout<<"Use u for car starting"<<endl;
    cout<<"Use U for car returning "<<endl;
    cout<<"Use '0' for light 0"<<endl;
    cout<<"Use '1' to stop light0"<<endl;
    cout<<"Use '2' for ambient light0"<<endl;
    cout<<"Use '3' for diffuse  light0"<<endl;
    cout<<"Use '4' for specular  light0"<<endl;

    cout<<"Use '5' for light 1"<<endl;
    cout<<"Use '6' to stop light1"<<endl;
    cout<<"Use '7' for ambient light1"<<endl;
    cout<<"Use '8' for diffuse  light1"<<endl;
    cout<<"Use '9' for specular  light1"<<endl;

    glutMainLoop();

    return EXIT_SUCCESS;
}
