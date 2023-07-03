/*! \mainpage Modélisation Géométrique d'un Maillage
 
 
 \section cont Contexte
 
 Dans \f$\mathbb{R}^3\f$, la construction d'une 1- ou d'une 2-variété est une opération courante en CAO.
 
 <b>Distinctions majeures d'un maillage à l'autre :</b> 
	- Nature de la maille élémentaire (détermination implicite de la structuration du maillage), 
	- Fermeture,
	- Genre.
  */


/**
 * \file main.cpp
 * \brief Programme principal
 * \author TB
 * \version 0.1
 * \date 09/10/18

 * Programme principal pour lancer l'interface du constructeur d'offset.
 */


//== CONSTANTS =================================================================

//#define OS_GENERIC
#define OS_WINDOWS
//#define OS_LINUX

#define PATH_DATA "dat/"
//#define FILENAME "bumpy_sphere.off"
#define FILENAME "rounded_cube1.stl"

//torus.off : 288v
//bumpy_sphere.off : 5724v
//bunny.off : 34835v
//sphere.off : 642v
//scanned_face.off : 110757v 


//== INCLUDES =================================================================

#include <iostream>
#include "Viewer.hpp"

#ifdef OS_GENERIC
    #include <ctime>
#endif
#ifdef OS_WINDOWS
    #include <windows.h>
#endif
#ifdef OS_LINUX
    #include <sys/time.h>
#endif


/**
 * \fn void Idle()
 * \brief  Fonction exécutée automatiquement par GLUT s'il n'y a pas d'événement en attente et donc si l'application est disponible pour la réalisation d'une tâche annexe. La tâche annexe en question ici est la fermeture de la console d'affichage si l'interface est fermée.
 */
void Idle(){

if(glutGet(GLUT_WINDOW_FORMAT_ID)==0) exit(0);

} 

/**
 * \fn string  get_file_name(string path_file_name, string _file_name)
 * \brief Format properly the filename of the mesh
 */
string  get_file_name(string  path_file_name, string _file_name){
	return path_file_name + _file_name;
}


/**
 * \fn int main(int argc, char **argv)
 * \brief Programme principal
 * \param argc Nombre de paramètres passé au programme. Ce nombre est toujours >1, car le premier est le nom de l'exécutable.
 * \param argv Tableau de chaînes de caractères contenant les paramètres passés au programme.
 * \return EXIT_SUCCESS - Arrêt normal du programme.
 */
int main(int argc, char **argv)
{


//VIEWER
//*************************

  glutInit(&argc, argv);

  //Create a viewer
  Viewer viewer("Mesh Viewer", 768, 512);

  cout<<endl<<"File"<<endl<<"****************"<<endl;

  //Read the CAO file
  if (argc>1){
    viewer.open_mesh(argv[1]);
  }
  else{

	  viewer.open_mesh(get_file_name(PATH_DATA, FILENAME).c_str());
  }
 
  glutIdleFunc(Idle); 


//GLUT LOOP
//*************************
  glutMainLoop();

  return EXIT_SUCCESS;
}