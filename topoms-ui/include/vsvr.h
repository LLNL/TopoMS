/**
 * @file    vsvr.h
 * @author  Thomas Lewiner <tomlew@mat.puc-rio.br>
 * @author  Math Dept, PUC-Rio
 * @version 2008.1
 * @date    07/03/2008
 *
 * @brief   Very Simple Volume Rendering
 */
/**
  Harsh Bhatia (hbhatia@llnl.gov) borrowed the online available code by
  Thomas Lewiner, and adapted to use for TopoMS
 */
// -----------------------------------------------------------------------------
#ifndef _VSVR_H_
#define _VSVR_H_

#if !defined(WIN32) || defined(__CYGWIN__)
#pragma interface
#endif // WIN32

#ifdef USE_GLEW
#include <GL/glew.h>
#endif

#ifdef __APPLE__
    #include <OpenGL/glu.h>  // openGL utilities
    #include <OpenGL/gl.h>   // openGL declarations
#else
 //   #include <GL/glu.h>  // openGL utilities
 // #include <GL/gl.h>   // openGL declarationd aniketh: glu.h already includes gl.h
#endif

#include <float.h>   // definition of FLT_EPSILON
#include <stdio.h>   // definition of printf
#include <stdlib.h>  // definition of NULL

//_____________________________________________________________________________
/** Very Simple Volume Rendering */
/** \class VSVR
  * \brief the Very Simple Volume Rendering containing a texture and transfer function.
  */
class VSVR
//-----------------------------------------------------------------------------
{
private:
    void glsl_log(int OBJ){

        int infologLength  = 0;
        glGetObjectParameterivARB(OBJ, GL_OBJECT_INFO_LOG_LENGTH_ARB, &infologLength);
        if (infologLength > 1) {

            int charsWritten  = 0;
            char * const log = new char[infologLength];
            glGetInfoLogARB(OBJ, infologLength, &charsWritten, log);
            printf(" VSVR() -- glsl_log: \n %s\n", log);
            delete log;
        }
    }


// Constructors
public :
  /**
   * Main and default constructor
   * \brief constructor
   * \param tex_ni width  of the 3D texture (must be a power of 2)
   * \param tex_nj depth  of the 3D texture (must be a power of 2)
   * \param tex_nk height of the 3D texture (must be a power of 2)
   * \param tf_size size of the transfer function
   */
    VSVR (const char* volShaderFilename = 0, const int tex_ni = -1, const int tex_nj = -1, const int tex_nk = -1, const int tf_size = -1 ) :
            _tex_extern (false), _tex_ni(tex_ni), _tex_nj(tex_nj), _tex_nk(tex_nk), _tex((float*)NULL),
            _tf_extern  (false), _tf_size(tf_size), _tf((float*)NULL), _rescale(true), _tex_glid(0) {


        //printf("VSVR::VSVR(%s)\n", volShaderFilename);

        const char* shader_source = (volShaderFilename == 0) ? get_shaderSource() : textFileRead(volShaderFilename);

        // create a fragment shader object
        GLhandleARB fragshader_object = glCreateShaderObjectARB(GL_FRAGMENT_SHADER_ARB);
        glShaderSourceARB(fragshader_object, 1, &shader_source, NULL);

        // compile the shader
        glCompileShaderARB(fragshader_object);

        // check compile status
        {
            int status;
            glGetShaderiv(fragshader_object, GL_COMPILE_STATUS, &status);
            if (status != GL_TRUE){
                glsl_log(fragshader_object);
            }
        }

        _sprog = glCreateProgramObjectARB();
        glAttachObjectARB(_sprog, fragshader_object);
        glLinkProgramARB(_sprog);

        // check link status
        {
            int status;
            glGetShaderiv(_sprog, GL_COMPILE_STATUS, &status);
            if (status != GL_TRUE){
                glsl_log(_sprog);
            }
        }
    }

  /** Destructor */
  ~VSVR() { tex_glunload() ;  tf_glunload() ;  tex_free() ;  tf_free() ; }


//-----------------------------------------------------------------------------
// Rendering
public :
  /**
   * loads the texture and renders
   * \param nslices number of slices cut into the texture
   * \param opacity factor to apply to the transfer function
   */
  bool gl_render    ( int nslices = 256 ) ;

  /**
   * redisplay with the current setting (possibly a different viewpoint)
   * \param nslices number of slices cut into the texture
   */
  bool gl_redisplay ( int nslices = 256 ) const ;

  /** accesses the texture name */
  int  tex_glid     () const { return _tex_glid ; }


protected :

    char *get_shaderSource() {

#ifdef WIN32
        return "uniform sampler3D voltex;\
                uniform sampler1D transf;\
                uniform float vmin;\
                uniform float vmax;\
                void main() {\
                float value = texture3D(voltex, gl_TexCoord[0].stp);\
                value = (value - vmin) / (vmax - vmin);\
                gl_FragColor = texture1D(transf, value);\
                }";
#else
        return "uniform sampler3D voltex;\
                uniform sampler1D transf;\
                uniform float vmin;\
                uniform float vmax;\
                void main() {\
                float value = texture3D(voltex, gl_TexCoord[0].stp).r;\
                value = (value - vmin) / (vmax - vmin);\
                gl_FragColor = texture1D(transf, value);\
                }";
#endif
    }

    char *textFileRead(const char *fn) {

        if(fn == NULL)
            return NULL;

        int count = 0;
        char *content = NULL;

        FILE *fp = fopen(fn,"rt");
        if(fp == NULL){
            printf(" Could not open file %s\n", fn);
            return NULL;
        }

        fseek(fp, 0, SEEK_END);
        count = ftell(fp);
        rewind(fp);

        if (count > 0) {
            content = (char *)malloc(sizeof(char) * (count+1));
            count = fread(content,sizeof(char),count,fp);
            content[count] = '\0';
        }
        fclose(fp);
        return content;
    }

  /** rescale the texels to match the transfer function size */
  void tex_rescale  () ;
  /** loads the 3D texture */
  bool tex_glload   () ;
  /**
   * loads the transfer function
   * \param opacity factor to apply to the transfer function
   */
  bool tf_glload    ()  ;

  /** unloads the 3D texture */
  void tex_glunload () ;
  /** unloads the transfer function */
  void tf_glunload  () ;

  /** sets the openGL attributes */
  void gl_set       () const ;
  /** unsets the openGL attributes */
  void gl_unset     () const ;

  /** sets the clipping planes (uses the 6 first) */
  void gl_clip      () const ;
  /** unsets the clipping planes */
  void gl_unclip    () const ;

/// -----------------------------------------------------------------------------
/// 3D texture accessors
public :
  /**  accesses the width  of the 3D texture */
  inline const int tex_ni() const { return _tex_ni ; }
  /**  accesses the depth  of the 3D texture */
  inline const int tex_nj() const { return _tex_nj ; }
  /**  accesses the height of the 3D texture */
  inline const int tex_nk() const { return _tex_nk ; }

  /**
   * changes the size of the 3D texture
   * \param tex_ni width  of the 3D texture (must be a power of 2)
   * \param tex_nj depth  of the 3D texture (must be a power of 2)
   * \param tex_nk height of the 3D texture (must be a power of 2)
   */
  inline void tex_set_resolution( const int tex_ni, const int tex_nj, const int tex_nk ) { _tex_ni = tex_ni ;  _tex_nj = tex_nj ;  _tex_nk = tex_nk ; }

  /*inline void tf_log(bool *tlog) {
      _tf_log = tlog;
  }*/

  /**
   * selects to use a 3D texture allocated from another class
   * \param tex is the pointer to the external 3D texture, allocated as a tex_ni*tex_nj*tex_nk vector running in i first. Its values will be rescaled.
   */
  inline void tex_set_extern ( float *tex )
  { tex_free() ;  _tex_extern = tex != NULL ;  _tex = tex ; }

  /**
   * selects to allocate the 3D texture
   */
  inline void tex_set_intern ()
  { if( _tex_extern ) _tex = NULL ;  _tex_extern = false ; }

  /** allocates the 3D texture */
  inline void tex_alloc  ()
  { tex_free() ;  int tex_size = _tex_ni*_tex_nj*_tex_nk ;  if( tex_size > 0 ) _tex = new float[tex_size] ;  _rescale = true ; }

  /** frees the 3D texture */
  inline void tex_free   ()
  { if( !_tex_extern ) delete [] _tex ; _tex = NULL ; }

  /**
   * accesses a specific voxel of the 3D texture
   * \param i abscisse of the voxel
   * \param j ordinate of the voxel
   * \param k height   of the voxel
   */
  inline const float tex_get ( const int i, const int j, const int k ) const
  { return _tex[ i + j*_tex_ni + k*_tex_nj*_tex_ni ] ; }

  /**
   * sets the value of a specific voxel of the 3D texture
   * \param val new value for the voxel
   * \param i abscisse of the voxel
   * \param j ordinate of the voxel
   * \param k height   of the voxel
   */
  inline void        tex_set ( const int i, const int j, const int k, const float val )
  { _tex[ i + j*_tex_ni + k*_tex_nj*_tex_ni ] = val ; }


/// -----------------------------------------------------------------------------
///  Transfer function (color map) accessors
public :

  /**  accesses the size of the transfer function */
  inline const int tf_size() const
  { return _tf_size ; }

  /**
   * changes the size of the transfer function
   * \param tf_size size of the transfer function
   */
  inline void tf_set_size ( const int tf_size )
  { if( _tf_size != tf_size ) _rescale = true ;  _tf_size = tf_size ; }

  /**
   * selects to use a transfer function  allocated from another class
   * \param tf is the pointer to the external data, allocated as a size_x*size_y*size_z vector running in x first
   */
  inline void tf_set_extern ( float *tf )
  { tf_free() ;  _tf_extern = tf != NULL ;  _tf = tf ; }

  /**
   * selects to allocate the transfer function
   */
  inline void tf_set_intern ()
  { if( _tf_extern ) _tf = NULL ;  _tf_extern = false ; }

  /** allocates the transfer function */
  inline void tf_alloc  ()
  { tf_free() ;  if( _tf_size > 0 ) _tf = new float[4*_tf_size] ; }

  /** frees the transfer function */
  inline void tf_free   ()
  { if( !_tf_extern ) delete [] _tf ; _tf = NULL ; }


  /**
   * accesses a specific element of the transfer function
   * \param i element index
   * \param r returned red   component of the color map
   * \param g returned green component of the color map
   * \param b returned blue  component of the color map
   * \param a returned transparency    of the color map
   */
  inline void tf_get ( const int i, float &r, float &g, float &b, float &a ) const
  { float *ptr = _tf + i ;  r = *ptr ;  ptr += tf_size() ;  g = *ptr ;  ptr += tf_size() ;  b = *ptr ;  ptr += tf_size() ;   a = *ptr ; }

  /**
   * sets a specific element of the transfer function
   * \param i element index
   * \param r red   component of the color map
   * \param g green component of the color map
   * \param b blue  component of the color map
   * \param a transparency    of the color map
   */
  inline void tf_set ( const int i, const float r, const float g, const float b, const float a )
  { float *ptr = _tf + i ;  *ptr = r ;  ptr += tf_size() ;  *ptr = g ;  ptr += tf_size() ;  *ptr = b ;  ptr += tf_size() ;  *ptr = a ; }


/// -----------------------------------------------------------------------------
/// Elements
protected :
  bool      _tex_extern ;  /**< selects wether to allocate the 3D texture or to use one allocated from another class */
  int       _tex_ni     ;  /**< width  of the 3D texture (must be a power of 2) */
  int       _tex_nj     ;  /**< depth  of the 3D texture (must be a power of 2) */
  int       _tex_nk     ;  /**< height of the 3D texture (must be a power of 2) */
  float    *_tex        ;  /**< the 3D texture : grid of float values */

  float _texmin, _texmax;
  //bool  *_tf_log;

  bool      _tf_extern  ;  /**< selects wether to allocate the transfer function or to use one allocated from another class */
  int       _tf_size    ;  /**< size of the transfer function */
  float    *_tf         ;  /**< the transfer function : colormap with 4 floats (rgba) per color*/

private :
    bool   _rescale    ;  /**< needs to rescale */
    GLuint _tf_glid    ;
    GLuint _tex_glid   ;  /**< openGL texture name */

    GLuint _sprog      ;  /**< the GLSL fragment shader */
};
/// _____________________________________________________________________________


#endif // _VSVR_H_
