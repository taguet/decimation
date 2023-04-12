
//=============================================================================
//
//  CLASS GlutViewer - IMPLEMENTATION
//
//=============================================================================

//== INCLUDES =================================================================

#include "gl.h"
#include "GlutViewer.hpp"


//== IMPLEMENTATION ========================================================== 


map<int, GlutViewer*>  GlutViewer::windows__;


//-----------------------------------------------------------------------------


GlutViewer::
GlutViewer(const char* _title, int _width, int _height) : width_(_width), height_(_height) {

  // create window
  glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE | GLUT_ALPHA);
  glutInitWindowSize(_width, _height);
  windowID_ = glutCreateWindow(_title);
  int mainID = windowID_;
  windows__[windowID_] = this;

  // register callbacks
  glutDisplayFunc(display__);
  glutKeyboardFunc(keyboard__);
  glutSpecialFunc(special__);
  glutMouseFunc(mouse__);
  glutMotionFunc(motion__);
  glutPassiveMotionFunc(passivemotion__);
  glutReshapeFunc(reshape__); 
  glutVisibilityFunc(visibility__);

  // create subwindow 
  /*
  windowID_ = glutCreateSubWindow(windowID_, 0, 0, 10, 20);
  windows__[windowID_] = this;
  glutDisplayFunc(subDisplay__);

  glutSetWindow(mainID);
  */

  // create menu
  n_draw_modes_ = 0;
  menuID_ = glutCreateMenu(processmenu__);
  glutAttachMenu(GLUT_RIGHT_BUTTON);

}
  

//-----------------------------------------------------------------------------


GlutViewer::
~GlutViewer()
{
  glutDestroyWindow(windowID_);
  glutDestroyMenu(menuID_);
}


//----------------------------------------------------------------------------


void
GlutViewer::clear_draw_modes()
{
  for (unsigned int i=0; i<n_draw_modes_; ++i)
    glutRemoveMenuItem(1);

  n_draw_modes_ = 0;
  draw_mode_names_.clear();
}


unsigned int
GlutViewer::add_draw_mode(const string& _s)
{
  // insert in popup menu
  glutAddMenuEntry(_s.c_str(), n_draw_modes_);

  ++n_draw_modes_;
  draw_mode_names_.push_back(_s);

  return n_draw_modes_-1;
}


void
GlutViewer::set_draw_mode(int _id)
{
  draw_mode_ = _id;
  glutPostRedisplay();
}


int
GlutViewer::get_draw_mode()
{
  return draw_mode_;
}


//-----------------------------------------------------------------------------


GlutViewer* GlutViewer::current_window() { 
  return windows__[glutGetWindow()]; 
}

void GlutViewer::display__(void) {
  current_window()->display();
}

void GlutViewer::subDisplay__(void) {
  current_window()->displaySub();
}

void GlutViewer::idle__(void) {
  current_window()->idle();
} 

void GlutViewer::keyboard__(unsigned char key, int x, int y) {
  current_window()->keyboard((int)key, x, y);
}

void GlutViewer::motion__(int x, int y) {
  current_window()->motion(x, y);
}

void GlutViewer::mouse__(int button, int state, int x, int y) {
  current_window()->mouse(button, state, x, y);
}

void GlutViewer::passivemotion__(int x, int y) {
  current_window()->passivemotion(x, y);
}

void GlutViewer::reshape__(int w, int h) {
  current_window()->reshape(w, h);
}

void GlutViewer::special__(int key, int x, int y) {
  current_window()->special(key, x, y);
}   

void GlutViewer::visibility__(int visible) {
  current_window()->visibility(visible);
}

void GlutViewer::processmenu__(int id) {
  current_window()->processmenu(id);
}


//-----------------------------------------------------------------------------


void GlutViewer::idle(void) {} 
void GlutViewer::motion(int x, int y) {}
void GlutViewer::mouse(int button, int state, int x, int y) {}
void GlutViewer::passivemotion(int x, int y) {}
void GlutViewer::visibility(int visible) {}
void GlutViewer::reshape(int w, int h) {}


void 
GlutViewer::processmenu(int i) 
{
  set_draw_mode(i); 
}


void GlutViewer::keyboard(int key, int x, int y) 
{
  switch (key)
  {
    case 27:  
    {
      exit(0); 
      break;
    }
	}
}


void GlutViewer::special(int key, int x, int y) 
{
  switch (key)
  {
    case GLUT_KEY_F12: 
    {
      if (!fullscreen_) 
      {
				bak_left_   = glutGet(GLUT_WINDOW_X);
				bak_top_    = glutGet(GLUT_WINDOW_Y);
				bak_width_  = glutGet(GLUT_WINDOW_WIDTH);
				bak_height_ = glutGet(GLUT_WINDOW_HEIGHT);
				glutFullScreen();
				fullscreen_ = true;
      }
      else
      {
				glutReshapeWindow(bak_width_, bak_height_);
				glutPositionWindow(bak_left_, bak_top_);
				fullscreen_ = false;
      }
      break;
    }
  }
} 


void GlutViewer::display(void) 
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  if (draw_mode_ < draw_mode_names_.size())
    draw(draw_mode_names_[draw_mode_]);
  else
    draw("");

  glutSwapBuffers();
}


void GlutViewer::displaySub(void)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    draw("");
    glutSwapBuffers();
}


//=============================================================================
