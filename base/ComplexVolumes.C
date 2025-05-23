#include "ComplexVolumes.h"

#include "Vector.h"
#include "Fields.h"
#include "VolumeGeometry.h"
#include "ObjParser.h"
#include <iostream>
#include <iomanip>
#include <cmath>


void makeObjVolume(std::string filename, ScalarField& density, Mesh& mesh, int gridRes, int bandwidth)
{
	ScalarGrid grid = ScalarGrid(new SGrid<float>);

	grid->setDefVal(0.0);

	// Parse obj
	ObjParser parser;
	mesh = createEmptyMesh();
	parser.load(filename, mesh);
	Vector gridDims = mesh->URC() - mesh->LLC();
	grid->init(gridRes, gridRes, gridRes, gridDims.X(), gridDims.Y(), gridDims.Z(), mesh->LLC() );
	
	// Ray march the level set of the mesh
	RayMarchLevelSet(mesh, grid, bandwidth);

	// Make grid into volume
	density = volumize(grid);
	
	float scale_amt = 1;
	density = scale(density, Vector(scale_amt,-scale_amt,scale_amt));

}

void makeTikeMyson(ScalarField& density, ColorField& color)
{

	float pi = 3.14159265359;


	// BODY
	Vector body_pos = -Vector();
	Vector body_dir = Vector(0,-1,0);
	float body_height = 0.6;
	float body_angle =  (35.0/360.0)*2*pi;

	ScalarField body = Cone(body_pos, body_dir, body_height, body_angle);	
	body = scale(body, Vector(1, 1, .6));

	// HEAD
	float head_rad = 0.2;
	Vector head_pos = -Vector(0, body_height + head_rad/2.0, 0);
	
	ScalarField head = Sphere(head_pos, head_rad);

	// TEETH
	Vector tooth_pos = head_pos + Vector(.05, .12, .2);
	
	ScalarField teeth = CsgBox(Vector(), .08, 8);
	teeth = scale(teeth, Vector(0.5,1,0.5));

	teeth = intersection(teeth, Plane(Vector(), Vector(0,-1,0)));
	teeth = translate(teeth, tooth_pos);
	teeth = Union(teeth, scale(teeth, Vector(-1,1,1)) );

	// EARS
	
	Vector ear_pos = -Vector(0.25,.95,0);
	float ear_rad = 0.18;

	ScalarField ear = Sphere(Vector(), ear_rad);
	ear = scale(ear, Vector(1,1,.5));
	ear = translate(ear, ear_pos);
	ear = Union(ear, scale(ear, Vector(-1,1,1)));
	
	head = Union(head, ear);
/*
	ScalarField measureStick = CsgBox(Vector(), 1, 4);
	measureStick = scale(measureStick, Vector(.6, .7, .4));
	measureStick = translate(measureStick, Vector(0.6,0,0));
	head = Union(head, measureStick);
*/
	// EYES
	
	float eye_size = 0.15;
	float eye_spacing = .075;
	float eye_forward = .2;
	float eye_height = head_pos.Y() - head_rad/4.0;

	ScalarField eye = scale(SteinerPatch(), Vector(eye_size, eye_size, eye_size));
	eye = Union( translate(eye, Vector(eye_spacing, eye_height, eye_forward)),  translate(eye, Vector(-eye_spacing, eye_height, eye_forward)) );

	// SHOULDERS
	
	Vector shoulder_pos = -Vector(0.37,.55,0);
	float shoulder_rad = 0.165;

	ScalarField shoulder = Sphere(Vector(), shoulder_rad);
	shoulder = translate(shoulder, shoulder_pos);
	shoulder = Union(shoulder, scale(shoulder, Vector(-1,1,1)));
	
	// ARMS
	float arm_length = body_height * (1.0/2.25);
	float arm_width = arm_length * .55;
	Vector upper_arm_pos = -Vector(.5, 2.0*body_height/3.0, .1);
	Vector upper_arm_axis = Vector(-.7,1,-upper_arm_pos.Z());
	
	ScalarField upper_left_arm = Ellipse(upper_arm_pos, upper_arm_axis, arm_length, arm_width);
	upper_arm_pos[0] *= -1;
	upper_arm_axis[0] *= -1;
	ScalarField upper_right_arm = Ellipse(upper_arm_pos, upper_arm_axis, arm_length, arm_width);

	Vector lower_arm_pos = -Vector(.5, 0, .1);
	Vector lower_arm_axis = Vector(.7,1,lower_arm_pos.Z());
	
	ScalarField lower_left_arm = Ellipse(lower_arm_pos, lower_arm_axis, arm_length, arm_width*2.0/3.0);
	lower_arm_pos[0] *= -1;
	lower_arm_axis[0] *= -1;
	ScalarField lower_right_arm = Ellipse(lower_arm_pos, lower_arm_axis, arm_length, arm_width*2.0/3.0);


	ScalarField left_arm = Union(upper_left_arm, lower_left_arm);
	ScalarField right_arm = Union(upper_right_arm, lower_right_arm);

	ScalarField arms = Union(left_arm, right_arm);

	// HANDS
	float handsize = 0.02;
	float handspace = .3;
	float handforward = -.1;
	float handheight = .2;

	ScalarField hand = Icosahedron();
	hand = scale(hand, Vector(handsize, handsize*1.3, handsize));
	hand = intersection(hand, Plane(Vector(), Vector(0,1,0)));
	hand = rotate(hand, 45.0/360.0 *2*pi, Vector(0,0,1));
	hand = Union( translate(hand, Vector(handspace, handheight, handforward)), translate(rotate(hand, pi, Vector(0,1,0)), Vector(-handspace, handheight, handforward)));
	hand = add(hand, constant(.03));
	arms = Union(arms, hand);

	// WRISTBANDS
	
	ScalarField wristband = Torus(Vector(.4, .13, -.125),Vector(-1,1,.2), .1, .03);
	wristband = Union(wristband, scale(wristband, Vector(-1,1,1)));
	
	// LEGS
	float leg_length = body_height * (1.0/2.25);
	float leg_width = leg_length * .35;
	Vector upper_leg_pos = Vector(.1, leg_length/2.0,0);
	Vector upper_leg_axis = Vector(0,1,0);
	
	ScalarField upper_left_leg = Ellipse(upper_leg_pos, upper_leg_axis, leg_length, leg_width);
	upper_leg_pos[0] *= -1;
	ScalarField upper_right_leg = Ellipse(upper_leg_pos, upper_leg_axis, leg_length, leg_width);

	Vector lower_leg_pos = Vector(.1, upper_leg_pos.Y()+leg_length, 0);
	Vector lower_leg_axis = Vector(0,1,0);
	
	ScalarField lower_left_leg = Ellipse(lower_leg_pos, lower_leg_axis, leg_length, leg_width*2.0/3.0);
	lower_leg_pos[0] *= -1;
	lower_leg_axis[0] *= -1;
	ScalarField lower_right_leg = Ellipse(lower_leg_pos, lower_leg_axis, leg_length, leg_width*2.0/3.0);


	ScalarField left_leg = Union(upper_left_leg, lower_left_leg);
	ScalarField right_leg = Union(upper_right_leg, lower_right_leg);

	ScalarField legs = Union(left_leg, right_leg);
	
	// FEET
	float foot_length =.17;
	float foot_width  = .07;
	Vector foot_pos = Vector(.1, 0,foot_length/2.0);
	Vector foot_axis = Vector(0,0,1);
	
	ScalarField left_foot = Ellipse(foot_pos, foot_axis, foot_length, foot_width);
	foot_pos[0] *= -1;
	ScalarField right_foot = Ellipse(foot_pos, foot_axis, foot_length, foot_width);

	ScalarField feet = Union(left_foot, right_foot);
	feet = scale(feet, Vector(1, .7, 1));
	feet = translate(feet, Vector(0, lower_leg_pos.Y()+leg_length, 0));

	// PANTS
	float pants_length = .1;
	float pants_roundness = 4;
	float pants_thickness = .01;
	Vector pants_pos(0,.1,0);
	
	ScalarField pants = CsgBox(Vector(), pants_length, pants_roundness);
	pants = scale(pants, Vector(.75,1,.5));
	pants = translate(pants,pants_pos);

	// remove pants from body
	body = cutout(body, pants);

	pants = Shell(pants, pants_thickness);
	pants = intersection(pants, Plane(pants_pos, Vector(0,1,0)));

		
	ScalarField belt = Torus(Vector(0,-.15,0),Vector(0,1,0), .15, .04);
	pants = Union(pants, belt);

	// MAIN DENSITIES
	float s = 0.3;
	ScalarField clothes = Union(wristband, pants);
	clothes = Union(clothes, teeth);
	clothes = scale(mask(clothes), Vector(s,s,s));

	ScalarField skin;
	skin = Union(body, head);
	skin = Union(skin, arms);
	skin = Union(skin, shoulder);
	skin = Union(skin, legs);
	skin = Union(skin, feet);
	skin = cutout(skin, eye);
	// sssssssss
	ScalarField skin_density = scale(mask(skin), Vector(s,s,s));
	eye = scale(multiply(mask(eye),2), Vector(s,s,s));

	// COLORS

	Color skin_color(1.0,0.5,.25,0);
	ColorField skin_colorField = constant(skin_color);
	skin_colorField = multiply(skin_colorField, skin_density);

	Color eye_color(1,0,0,0);
	ColorField eye_colorField = multiply(constant(eye_color), constant(1.5));
	eye_colorField = multiply(eye_colorField, eye);


	Color clothes_color(1,1,1,0);
	ColorField clothes_colorField = multiply(constant(clothes_color), constant(1.5));
	clothes_colorField = multiply(clothes_colorField, clothes);


	color = add(skin_colorField, eye_colorField);
	color = add(color, clothes_colorField);

	density = skin_density;
	density = Union(skin_density, eye);
	density = Union(density, clothes);
	density = multiply(density, 4);
	
}
