                                                                                                                    // POV-Ray 3.7 Scene File " ... .pov"
// author:  ...
// date:    ...
//------------------------------------------------------------------------
#version 3.7;
global_settings{ assumed_gamma 1.0 }
#default{ finish{ ambient 0.1 diffuse 0.9 }} 
//------------------------------------------------------------------------
#include "colors.inc"
#include "textures.inc"
#include "glass.inc"
#include "metals.inc"
#include "golds.inc"
#include "stones.inc"
#include "woods.inc"
#include "shapes.inc"
#include "shapes2.inc"
#include "functions.inc"
#include "math.inc"
#include "transforms.inc"
//------------------------------------------------------------------------
#declare Camera_0 = camera {/*ultra_wide_angle*/ angle 15      // front view
                            location  <0.0 , 1.0 ,-400.0>
                            right     x*image_width/image_height
                            look_at   <0.0 , 1.0 , 0.0>}
#declare Camera_1 = camera {/*ultra_wide_angle*/ angle 30   // diagonal view
                            location  <10.0 , 5.0 ,-10.0>
                            right     x*image_width/image_height
                            look_at   <0.0 , 1 , 0.0>}
#declare Camera_2 = camera {/*ultra_wide_angle*/ angle 90  //right side view
                            location  <3.0 , 1.0 , 0.0>
                            right     x*image_width/image_height
                            look_at   <0.0 , 1.0 , 0.0>}
#declare Camera_3 = camera {/*ultra_wide_angle*/ angle 90        // top view
                            location  <0.0 , 3.0 ,-0.001>
                            right     x*image_width/image_height
                            look_at   <0.0 , 1.0 , 0.0>}   
#declare Camera_4 = camera {/*ultra_wide_angle*/ angle 90        // top view
                            location  <4, 12 ,6>
                            right     x*image_width/image_height
                            look_at   <4, 1.0 , 6>}       
#declare Camera_5 = camera {/*ultra_wide_angle*/ angle 90        // top view for longer actuators
                            location  <6, 13 ,8>
                            right     x*image_width/image_height
                            look_at   <6, 1.0 , 8>}                                                   
camera{Camera_5}
//------------------------------------------------------------------------
// sun -------------------------------------------------------------------
light_source{<0,2500,0> color White}
// sky -------------------------------------------------------------------
sky_sphere{ pigment{ gradient <0,1,0>
                     color_map{ [0   color rgb<1,1,1>         ]//White
                                [0.4 color rgb<0.14,0.14,0.56>]//~Navy
                                [0.6 color rgb<0.14,0.14,0.56>]//~Navy
                                [1.0 color rgb<1,1,1>         ]//White
                              }
                     scale 2 }
           } // end of sky_sphere 
//------------------------------------------------------------------------

//------------------------------ the Axes --------------------------------
//------------------------------------------------------------------------
#macro Axis_( AxisLen, Dark_Texture,Light_Texture) 
 union{
    cylinder { <0,0,0>,<AxisLen,0,0>,0.1
               texture{checker texture{Dark_Texture } 
                               texture{Light_Texture}
                       translate<0.1,0,0.1>}
             }
    cone{<AxisLen,0,0>,0.4,<AxisLen+0.7,0,0>,0
          texture{Dark_Texture}
         }
     } // end of union                   
#end // of macro "Axis()"
//------------------------------------------------------------------------
#macro AxisXYZ( AxisLenX, AxisLenY, AxisLenZ, Tex_Dark, Tex_Light)
//--------------------- drawing of 3 Axes --------------------------------
#local Font_Name = "cyrvetic.ttf" 
union{
#if (AxisLenX != 0)
 object { Axis_(AxisLenX, Tex_Dark, Tex_Light)   rotate< 0,0,00> no_shadow}// x-Axis
 text   { ttf Font_Name,  "x",  0.15,  0  texture{Tex_Dark} 
          rotate<90,0,0> scale 1.5 translate <AxisLenX+1,0.4,-0.50> no_shadow}
#end // of #if 
#if (AxisLenY != 0)
 object { Axis_(AxisLenY, Tex_Dark, Tex_Light)   rotate< 0,0,  90> no_shadow}// y-Axis
 text   { ttf Font_Name,  "y",  0.15,  0  texture{Tex_Dark}    
          rotate<0,0,0> scale 0.75 translate <-0.65,AxisLenY+0.50,-0.10>  rotate<0,-45,0> no_shadow}
#end // of #if 
#if (AxisLenZ != 0)
 object { Axis_(AxisLenZ, Tex_Dark, Tex_Light)   rotate<0,-90, 0> no_shadow}// z-Axis
 text   { ttf Font_Name,  "z",  0.15,  0  texture{Tex_Dark} 
          rotate<90,0,0> scale 1.5 translate <-1,0.2,AxisLenZ+1> no_shadow}
#end // of #if 
} // end of union
#end// of macro "AxisXYZ( ... )"
//------------------------------------------------------------------------

#declare Texture_A_Dark  = texture {
                               pigment{ color Black}
                               finish { phong 1}
                             }
#declare Texture_A_Light = texture { 
                               pigment{ color Red}
                               finish { phong 1}
                             }

object{ AxisXYZ( 10, 0, 10, Texture_A_Dark, Texture_A_Light)}
//-------------------------------------------------- end of coordinate axes
 
 
 

// ground -----------------------------------------------------------------
//---------------------------------<<< settings of squared plane dimensions
#declare RasterScale = 1.0;
#declare RasterHalfLine  = 0.025*2;  
#declare RasterHalfLineZ = 0.025*2; 
//-------------------------------------------------------------------------
#macro Raster(RScale, HLine) 
       pigment{ gradient x scale RScale
                color_map{[0.000   color rgbt<1,1,1,0>*1.0]
                          [0+HLine color rgbt<1,1,1,0>*1.0]
                          [0+HLine color rgbt<1,1,1,1>]
                          [1-HLine color rgbt<1,1,1,1>]
                          [1-HLine color rgbt<1,1,1,0>*1.0]
                          [1.000   color rgbt<1,1,1,0>*1.0]} }
 #end// of Raster(RScale, HLine)-macro    
//-------------------------------------------------------------------------
    

plane { <0,1,0>, 0    // plane with layered textures
        texture { pigment{color MediumGoldenrod} }
        texture { Raster(RasterScale,RasterHalfLine ) rotate<0,0,0> }
        texture { Raster(RasterScale,RasterHalfLineZ) rotate<0,90,0>}
        rotate<0,0,0>
        translate<0,-5,0>
      }
//------------------------------------------------ end of squared plane XZ

//--------------------------------------------------------------------------
//---------------------------- objects in scene ----------------------------
//--------------------------------------------------------------------------

