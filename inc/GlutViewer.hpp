
//=============================================================================
//
//  CLASS GlutViewer
//
//=============================================================================


#ifndef GLUTVIEWER_HH
#define GLUTVIEWER_HH


//== INCLUDES =================================================================

#include "gl.h"
#include <map>
#include <vector>
#include <string>


//== CLASS DEFINITION =========================================================

	      

/** \class GlutViewer GlutViewer.hpp
    Simple Glut viewer. 
    Based on C++ glut interface of George Stetten and Korin Crawford.
**/

class GlutViewer
{
public:
   
  GlutViewer(const char* _title, int _width, int _height);
  virtual ~GlutViewer();



protected:

  virtual void draw(const std::string& _drawmode) = 0;
  void clear_draw_modes();
  unsigned int add_draw_mode(const std::string& _s);

  void set_draw_mode(int _id);  
  int get_draw_mode();

  void toggle_idle(bool _b) { glutIdleFunc( _b ? idle__ : NULL );	}

  virtual void display(void);
  virtual void displaySub(void);
  virtual void idle(void); 
  virtual void keyboard(int key, int x, int y);
  virtual void special(int key, int x, int y);
  virtual void motion(int x, int y);
  virtual void mouse(int button, int state, int x, int y);
  virtual void passivemotion(int x, int y);
  virtual void reshape(int w, int h); 
  virtual void visibility(int visible);
  virtual void processmenu(int i);

  int  width_, height_;


private:

  static void display__(void);
  static void subDisplay__(void);
  static void idle__(void); 
  static void keyboard__(unsigned char key, int x, int y);
  static void motion__(int x, int y);
  static void mouse__(int button, int state, int x, int y);
  static void passivemotion__(int x, int y);
  static void reshape__(int w, int h); 
  static void special__(int key, int x, int y);   
  static void visibility__(int visible);
  static void processmenu__(int i);

  static std::map<int, GlutViewer*>  windows__;
  static GlutViewer* current_window();

  

private:

  int  windowID_, menuID_; 

  bool fullscreen_;
  int  bak_left_, bak_top_, bak_width_, bak_height_;

  unsigned int              draw_mode_;
  unsigned int              n_draw_modes_;
  std::vector<std::string>  draw_mode_names_;
};


//=============================================================================
#endif // ACG_GLUTVIEWER_HH defined
//=============================================================================

