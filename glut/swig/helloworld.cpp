#include <GL/glut.h>


void renderScene(void) {
	glClear(GL_COLOR_BUFFER_BIT);
	glBegin(GL_TRIANGLES);
		glVertex3f(-0.5,-0.5,0.0);
		glVertex3f(0.5,0.0,0.0);
		glVertex3f(0.0,0.5,0.0);
	glEnd();
	glFlush();
}

void dome() {
  int argc=1;
  char *n = "inSight";
  char **argv=&n;
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_SINGLE | GLUT_RGBA);
	glutInitWindowPosition(100,100);
	glutInitWindowSize(320,320);
	glutCreateWindow("3D Tech- GLUT Tutorial");
	glutDisplayFunc(renderScene);
	glutMainLoop();
}

int main(int argc, char **argv) {
  dome();
	return 0;
}

