#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<windows.h>
#include<glut.h>
#include<vector>
#include<iostream>
#include<fstream>
#include<string>
using namespace std;

#define pi (2*acos(0.0))
#define cosA 0.99863
#define sinA 0.05234
#define sizeOfObject 20
#define FOV_Y 80
#define WINDOW_SIZE 500

double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
double angle;
double objectSize;

int level_of_recursion;
int num_of_pixel;
int num_of_object;
int num_of_lights;

#include "bitmap_image.hpp"
#include "FILE2.h"

struct point
{
	double x,y,z;
};

struct point cameraPos;
struct point u,r,l;

void drawAxes(){
	if(drawaxes==1)
	{
		glBegin(GL_LINES);{
			glColor3f(1.0, 0.0, 0.0);
			glVertex3f( 500,0,0);
			glVertex3f(-500,0,0);

			glColor3f(0.0, 1.0, 0.0);
			glVertex3f(0,-500,0);
			glVertex3f(0, 500,0);

			glColor3f(0.0, 0.0, 1.0);
			glVertex3f(0,0, 500);
			glVertex3f(0,0,-500);
		}glEnd();
	}
}

void drawlight(Vector3 pos){
	glColor3f(255/255.0, 20/255.0, 147/255.0);
	glPushMatrix();{
		glTranslatef(pos.x,pos.y,pos.z);
		glutSolidSphere(1,10, 10);
	}
	glPopMatrix();

}

bitmap_image image;
vector<vector<double> > z_buffer;
void capture(){
	cout << "in capture\n";
	image.setwidth_height(num_of_pixel, num_of_pixel);
	for(int i = 0; i < num_of_pixel; i++){
		for(int j = 0; j < num_of_pixel; j++){
			image.set_pixel(i,j,0,0,0);
		}
	}

	double plane_distance = (WINDOW_SIZE/2) / tan(pi*FOV_Y/360.0);
	Vector3 topleft;
	topleft.x = cameraPos.x + l.x * plane_distance - r.x * WINDOW_SIZE / 2 + u.x * WINDOW_SIZE / 2;
	topleft.y = cameraPos.y + l.y * plane_distance - r.y * WINDOW_SIZE / 2 + u.y * WINDOW_SIZE / 2;
	topleft.z = cameraPos.z + l.z * plane_distance - r.z * WINDOW_SIZE / 2 + u.z * WINDOW_SIZE / 2;

	double du = WINDOW_SIZE / float(num_of_pixel);
	double dv = du;

	cout << topleft.x << " " << topleft.y << " " << topleft.z << " " << plane_distance << " " <<  du << endl;
	for (int i = 0; i < num_of_pixel; ++i)
	{
		for (int j = 0; j < num_of_pixel; ++j)
		{
			Vector3 corner, eye;
			corner.x = topleft.x + j * r.x * du - i * u.x * dv - cameraPos.x;
			corner.y = topleft.y + j * r.y * du - i * u.y * dv - cameraPos.y;
			corner.z = topleft.z + j * r.z * du - i * u.z * dv - cameraPos.z;
			double norm = pow(corner.x * corner.x + corner.y * corner.y + corner.z * corner.z, 0.5);
			corner.x /= norm;  corner.y /= norm; corner.z /= norm;
			eye.x = cameraPos.x; eye.y = cameraPos.y; eye.z = cameraPos.z;

			Ray ray(eye, corner);

			double nearest = -1;
			double dummyColorAt[3];
			double t_min = INT_MAX;
			for (unsigned int k = 0; k < objects.size(); ++k)
			{
				double t = objects[k]->intersect(&ray, dummyColorAt, 0);
				if(t < 0) continue;
				if (t < t_min){
					t_min = t;
					nearest = k;
				}
			}

			if (nearest != -1)
			{
				double t = objects[nearest]->intersect(&ray, dummyColorAt, 1);
				image.set_pixel(j, i, 255 * dummyColorAt[0], 255 * dummyColorAt[1], 255 * dummyColorAt[2]);
			}
		}
	}

    cout << "capture done\n" ;
	image.save_image("capture.bmp");

}

void keyboardListener(unsigned char key, int x,int y){
	struct point temp;
	double ampt ;
	switch(key){
		case '1':
			temp.x = r.x ; temp.y = r.y ; temp.z = r.z ;
			r.x = r.x*cosA + l.x*sinA ;
			r.y = r.y*cosA + l.y*sinA ;
			r.z = r.z*cosA + l.z*sinA ;
			ampt = pow(r.x*r.x+r.y*r.y+r.z*r.z, 0.5);
			r.x = r.x/ampt ;
			r.y = r.y/ampt ;
			r.z = r.z/ampt ;

			l.x = l.x*cosA - temp.x*sinA ;
			l.y = l.y*cosA - temp.y*sinA ;
			l.z = l.z*cosA - temp.z*sinA ;
			ampt = pow(l.x*l.x+l.y*l.y+l.z*l.z, 0.5);
			l.x = l.x/ampt ;
			l.y = l.y/ampt ;
			l.z = l.z/ampt ;
			break;
		case '2':
			temp.x = r.x ; temp.y = r.y ; temp.z = r.z ;
			r.x = r.x*cosA - l.x*sinA ;
			r.y = r.y*cosA - l.y*sinA ;
			r.z = r.z*cosA - l.z*sinA ;
			ampt = pow(r.x*r.x+r.y*r.y+r.z*r.z, 0.5);
			r.x = r.x/ampt ;
			r.y = r.y/ampt ;
			r.z = r.z/ampt ;

			l.x = l.x*cosA + temp.x*sinA ;
			l.y = l.y*cosA + temp.y*sinA ;
			l.z = l.z*cosA + temp.z*sinA ;
			ampt = pow(l.x*l.x+l.y*l.y+l.z*l.z, 0.5);
			l.x = l.x/ampt ;
			l.y = l.y/ampt ;
			l.z = l.z/ampt ;
			break;
		case '4':
			temp.x = l.x ; temp.y = l.y ; temp.z = l.z ;
			l.x = l.x*cosA - u.x*sinA ;
			l.y = l.y*cosA - u.y*sinA ;
			l.z = l.z*cosA - u.z*sinA ;
			ampt = pow(l.x*l.x+l.y*l.y+l.z*l.z, 0.5);
			l.x = l.x/ampt ;
			l.y = l.y/ampt ;
			l.z = l.z/ampt ;

			u.x = u.x*cosA + temp.x*sinA ;
			u.y = u.y*cosA + temp.y*sinA ;
			u.z = u.z*cosA + temp.z*sinA ;
			ampt = pow(u.x*u.x+u.y*u.y+u.z*u.z, 0.5);
			u.x = u.x/ampt ;
			u.y = u.y/ampt ;
			u.z = u.z/ampt ;
			break;
		case '3':
			temp.x = l.x ; temp.y = l.y ; temp.z = l.z ;
			l.x = l.x*cosA + u.x*sinA ;
			l.y = l.y*cosA + u.y*sinA ;
			l.z = l.z*cosA + u.z*sinA ;
			ampt = pow(l.x*l.x+l.y*l.y+l.z*l.z, 0.5);
			l.x = l.x/ampt ;
			l.y = l.y/ampt ;
			l.z = l.z/ampt ;

			u.x = u.x*cosA - temp.x*sinA ;
			u.y = u.y*cosA - temp.y*sinA ;
			u.z = u.z*cosA - temp.z*sinA ;
			ampt = pow(u.x*u.x+u.y*u.y+u.z*u.z, 0.5);
			u.x = u.x/ampt ;
			u.y = u.y/ampt ;
			u.z = u.z/ampt ;
			break;
		case '5':
			temp.x = u.x ; temp.y = u.y ; temp.z = u.z ;
			u.x = u.x*cosA + r.x*sinA ;
			u.y = u.y*cosA + r.y*sinA ;
			u.z = u.z*cosA + r.z*sinA ;
			ampt = pow(u.x*u.x+u.y*u.y+u.z*u.z, 0.5);
			u.x = u.x/ampt ;
			u.y = u.y/ampt ;
			u.z = u.z/ampt ;

			r.x = r.x*cosA - temp.x*sinA ;
			r.y = r.y*cosA - temp.y*sinA ;
			r.z = r.z*cosA - temp.z*sinA ;
			ampt = pow(r.x*r.x+r.y*r.y+r.z*r.z, 0.5);
			r.x = r.x/ampt ;
			r.y = r.y/ampt ;
			r.z = r.z/ampt ;
			break;
		case '6':
			temp.x = u.x ; temp.y = u.y ; temp.z = u.z ;
			u.x = u.x*cosA - r.x*sinA ;
			u.y = u.y*cosA - r.y*sinA ;
			u.z = u.z*cosA - r.z*sinA ;
			ampt = pow(u.x*u.x+u.y*u.y+u.z*u.z, 0.5);
			u.x = u.x/ampt ;
			u.y = u.y/ampt ;
			u.z = u.z/ampt ;

			r.x = r.x*cosA + temp.x*sinA ;
			r.y = r.y*cosA + temp.y*sinA ;
			r.z = r.z*cosA + temp.z*sinA ;
			ampt = pow(r.x*r.x+r.y*r.y+r.z*r.z, 0.5);
			r.x = r.x/ampt ;
			r.y = r.y/ampt ;
			r.z = r.z/ampt ;
			break;
		case '0':
			capture();
			break;

		default:
			break;
	}
}

void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:	
			cameraPos.x -= 2 * l.x ;
			cameraPos.y -= 2 * l.y ;
			cameraPos.z -= 2 * l.z ;
			break;
		case GLUT_KEY_UP:		
			cameraPos.x += 2 * l.x ;
			cameraPos.y += 2 * l.y ;
			cameraPos.z += 2 * l.z ;
			break;

		case GLUT_KEY_RIGHT:
			cameraPos.x += 2 * r.x ;
			cameraPos.y += 2 * r.y ;
			cameraPos.z += 2 * r.z ;
			break;
		case GLUT_KEY_LEFT:
			cameraPos.x -= 2 * r.x ;
			cameraPos.y -= 2 * r.y ;
			cameraPos.z -= 2 * r.z ;
			break;

		case GLUT_KEY_PAGE_UP:
			cameraPos.x += 2 * u.x ;
			cameraPos.y += 2 * u.y ;
			cameraPos.z += 2 * u.z ;
			break;
		case GLUT_KEY_PAGE_DOWN:
			cameraPos.x -= 2 * u.x ;
			cameraPos.y -= 2 * u.y ;
			cameraPos.z -= 2 * u.z ;
			break;

		case GLUT_KEY_INSERT:
			break;
		case GLUT_KEY_HOME:
			if(objectSize > 0) objectSize -= 1 ;
			break;
		case GLUT_KEY_END:
			if(objectSize < 20) objectSize += 1 ;
			break;

		default:
			break;
	}
}

void mouseListener(int button, int state, int x, int y){	
	switch(button){
		case GLUT_LEFT_BUTTON:
			if(state == GLUT_DOWN){	
				drawaxes=1-drawaxes;
			}
			break;
		case GLUT_RIGHT_BUTTON:
			break;
		case GLUT_MIDDLE_BUTTON:
			break;
		default:
			break;
	}
}

void display(){
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);

	glLoadIdentity();
	gluLookAt(cameraPos.x,cameraPos.y,cameraPos.z,	cameraPos.x+10*l.x,cameraPos.y+10*l.y,cameraPos.z+10*l.z,	u.x,u.y,u.z);
	glMatrixMode(GL_MODELVIEW);

	drawAxes();
	for (unsigned int i=0; i < objects.size(); i++) {
		objects[i]->draw();
	}
	for (unsigned int i=0; i < lights.size(); i++) {
		drawlight(lights[i]);
	}
	glutSwapBuffers();
}


void animate(){
	angle+=0.05;
	glutPostRedisplay();
}

void init(){
	drawgrid=1;
	drawaxes=1;
	cameraHeight=150.0;
	cameraAngle=1.0;
	angle=0;
	objectSize = sizeOfObject/2;
	{
		cameraPos.x = 0;
		cameraPos.y = -130;
		cameraPos.z = 65;
	}
	{
		u.x = 0 ; u.y = 0 ; u.z = 1 ;
		r.x = 1 ; r.y = 0 ; r.z = 0 ;
		l.x = 0 ; l.y = 1 ; l.z = 0 ;

	}

	glClearColor(0,0,0,0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(FOV_Y,	1,	1,	1000.0);
}


void loadTestData(){
	Object *temp;
	Vector3 Center(0,0,10);
	double Radius = 10;
	temp=new Sphere(Center, Radius); 
	temp->setColor(1,0,0);
	temp->setCoEfficients(0.4,0.2,0.2,0.2);
	temp->setShine(1);
	objects.push_back(temp);

	/** temp=new Floor(1000, 20);
	temp->setCoEfficients(0.4,0.2,0.2,0.2);
	temp->setShine(1);
	objects.push_back(temp); */


	Vector3 light1(-50,50,50);
	lights.push_back(light1);

	num_of_pixel = 768;


}
void loadActualData(){
	ifstream infile;
	infile.open("scene.txt");
	infile >> level_of_recursion >> num_of_pixel >> num_of_object;
	for (int i = 0; i < num_of_object; ++i)
	{
		string line;
		infile >> line;
		if (line == "sphere")
		{
			double x,y,z, rad, red, green, blue, amb, diff, spec, reflec, shine;
			infile >> x >> y >> z >> rad >> red >> green >> blue >> amb >> diff >> spec >> reflec >> shine;
			Object *temp;
			Vector3 Center(x,y,z);
			temp=new Sphere(Center, rad);
			temp->setColor(red, green, blue);
			temp->setCoEfficients(amb, diff, spec, reflec);
			temp->setShine(shine);
			objects.push_back(temp);
		}
		if (line == "triangle")
		{
			Vector3 vect1, vect2, vect3;
			double red, green, blue, amb, diff, spec, reflec, shine;
			infile >> vect1.x >> vect1.y >> vect1.z >> vect2.x >> vect2.y >> vect2.z >> vect3.x >> vect3.y >> vect3.z;
			infile >> red >> green >> blue >> amb >> diff >> spec >> reflec >> shine;
			Object *temp;
			temp=new Triangle(vect1, vect2, vect3);
			temp->setColor(red, green, blue);
			temp->setCoEfficients(amb, diff, spec, reflec);
			temp->setShine(shine);
			objects.push_back(temp);

		}
		if(line == "general")
		{
			double A, B, C, D, E, F, G, H, I, J, length, width, height;
			double red, green, blue, amb, diff, spec, reflec, shine;
			Vector3 cube_ref_point;
			infile >> A >> B >> C >> D >> E >> F >> G >> H >> I >> J;
			infile >> cube_ref_point.x >> cube_ref_point.y >> cube_ref_point.z >> length >> width >> height;
			infile >> red >> green >> blue >> amb >> diff >> spec >> reflec >> shine;
			Object *temp;
			temp=new General(A, B, C, D, E, F, G, H, I, J);
			temp->setColor(red, green, blue);
			temp->setCoEfficients(amb, diff, spec, reflec);
			temp->setShine(shine);
			temp->setQuadParam(length, width, height);
			temp->setRefPoint(cube_ref_point);
			objects.push_back(temp);

		}
	}

	infile >> num_of_lights;
	for (int i = 0; i < num_of_lights; ++i)
	{
		Vector3 light;
		infile >> light.x >> light.y >> light.z;
		lights.push_back(light);
	}
	infile.close();

	Object *temp;
	temp=new Floor(1000, 20);
	temp->setCoEfficients(0.4,0.2,0.2,0.2);
	temp->setShine(1);
	objects.push_back(temp);

}

int main(int argc, char **argv){
	//loadTestData();
	loadActualData();

	glutInit(&argc,argv);
	glutInitWindowSize(WINDOW_SIZE, WINDOW_SIZE);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);

	glutCreateWindow("My OpenGL Program");

	init();

	glEnable(GL_DEPTH_TEST);

	glutDisplayFunc(display);
	glutIdleFunc(animate);

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();	

	return 0;
}
