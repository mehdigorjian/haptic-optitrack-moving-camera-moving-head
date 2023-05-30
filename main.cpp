#include <HDU/hduMath.h>
#include <HDU/hduMatrix.h>
#include <QHHeadersGLUT.h>  //Include all necessary headers

#include <cstdio>
///
#include <cmath>
#include <eigen3/Eigen/Geometry>

#include "cameras.h"
///
#include "OptiTrack.h"
///
#include <eigen3/Eigen/Geometry>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <thread>  // multithreading (using libs=-lpthread)
#include <vector>

#include "Model.h"
#include "Object.h"
#include "coordinateTransform.h"
#include "omp.h"  // parallel for loop (using libs=-lgomp and CXXFLAGS=-fopenmp)
//
#include "reactphysics3d/reactphysics3d.h"
Eigen::Vector3d hcurr_temp;
Eigen::Vector3d ocurr_temp;

std::vector<Eigen::Vector3d> optitrack_ptsA;
std::vector<Eigen::Vector3d> haptic_ptsB;
bool matrix_flag = false;

std::shared_ptr<coordinateTransform> t;  // = std::make_shared<coordinateTransform>();
std::shared_ptr<OptiTrack> opti = std::make_shared<OptiTrack>();

TriMesh* obj1 = nullptr;
TriMesh* obj2 = nullptr;
TriMesh* obj2_shadow = nullptr;
TriMesh* obj_coll_sp = nullptr;

Text* notCallibrate = nullptr;
Text* calibInstruction = nullptr;
////////////////////////////// physics engine
reactphysics3d::Vector3 updated_opti_obj_pos(0.0, 0.0, 0.0);
reactphysics3d::Vector3 updated_hpti_obj_pos(0.0, 0.0, 0.0);

reactphysics3d::Quaternion opti_obj_q;
reactphysics3d::Quaternion hpti_obj_q;
bool isOverlap;
bool isInside = false;
static float f_decay_count = 1.0f;

reactphysics3d::CollisionBody* body1_check;
reactphysics3d::CollisionBody* body2_check;
//////////////////////////////
hduVector3Dd eye(0., 0., 0.);

class DataTransportClass  // This class carried data into the ServoLoop thread
{
   public:
    TriMesh* c1;
    TriMesh* c2;
    TriMesh* c2_shadow;
    TriMesh* coll_sp;
};

double chargeRadius = 15.0;  // This variable defines the radius around the charge when the inverse square law changes to a spring force law.
double multiplierFactor = 40.0;
hduMatrix WorldToDevice;  // This matrix contains the World Space to DeviceSpace Transformation
hduVector3Dd forceVec;    // This variable contains the force vector.

int object_tracker = 1, screen_tracker = 3, head_tracker = 2;

// functions signatures
void physics(int argc, char** argv);
void glut_main(int, char**);
void button1DownCallback(unsigned int ShapeID);
void button1UpCallback(unsigned int ShapeID);
void touchCallback(unsigned int ShapeID);
void graphicsCallback(void);
void glutMenuFunction(int);
void HLCALLBACK computeForceCB(HDdouble force[3], HLcache* cache, void* userdata);                    // Servo loop callback
void HLCALLBACK startEffectCB(HLcache* cache, void* userdata);                                        // Servo Loop callback
void HLCALLBACK stopEffectCB(HLcache* cache, void* userdata);                                         // Servo Loop callback
hduVector3Dd forceField(hduVector3Dd Pos1, hduVector3Dd Pos2, HDdouble Multiplier, HLdouble Radius);  // This function computer the force beween the Model and the particle based on the positions

int main(int argc, char** argv) {
    std::thread t1(&OptiTrack::run, opti, argc, argv);
    std::thread t2(glut_main, argc, argv);
    std::thread t3(physics, argc, argv);
    t1.join();
    t2.join();
    t3.join();
    return 0;
}

void physics(int argc, char** argv) {
    reactphysics3d::PhysicsCommon physicsCommon;
    reactphysics3d::PhysicsWorld* world = physicsCommon.createPhysicsWorld();

    // const reactphysics3d::Vector3 halfExtents(5.0, 5.0, 5.0);

    reactphysics3d::Transform transform_opti_obj(updated_opti_obj_pos, opti_obj_q);
    reactphysics3d::Transform transform_hpti_obj(updated_hpti_obj_pos, hpti_obj_q);

    reactphysics3d::CollisionBody* body1 = world->createCollisionBody(transform_opti_obj);
    float radius1 = 5.0f;
    reactphysics3d::SphereShape* sphereShape1 = physicsCommon.createSphereShape(radius1);
    // reactphysics3d::BoxShape* sphereShape1 = physicsCommon.createBoxShape(halfExtents);
    reactphysics3d::Collider* collider1 = body1->addCollider(sphereShape1, transform_opti_obj);

    reactphysics3d::CollisionBody* body2 = world->createCollisionBody(transform_hpti_obj);
    float radius2 = 5.0f;
    reactphysics3d::SphereShape* sphereShape2 = physicsCommon.createSphereShape(radius2);
    // reactphysics3d::BoxShape* sphereShape2 = physicsCommon.createBoxShape(halfExtents);
    reactphysics3d::Collider* collider2 = body2->addCollider(sphereShape2, transform_hpti_obj);

    const reactphysics3d::decimal timeStep = 1.0f / 1000.0f;
    while (true) {
        world->update(timeStep);

        reactphysics3d::Transform newTransform_opti(updated_opti_obj_pos, opti_obj_q);
        reactphysics3d::Transform newTransform_hpti(updated_hpti_obj_pos, hpti_obj_q);
        body1->setTransform(newTransform_opti);
        body2->setTransform(newTransform_hpti);

        body1_check = body1;
        body2_check = body2;

        isOverlap = world->testOverlap(body1, body2);
        // if (isOverlap){
        //     std::cout << "-------------------------------------Overlap!" << std::endl;

        // }
    }
}

void glut_main(int argc, char** argv) {
    QHGLUT* DisplayObject = new QHGLUT(argc, argv);  // create a display window
    DisplayObject->setName("screen");
    DeviceSpace* deviceSpace = new DeviceSpace;  // Find a Phantom device named "Default PHANToM"
    DisplayObject->tell(deviceSpace);            // tell Quickhaptics that Omni exists
    DisplayObject->setBackgroundColor(0.0, 0.0, 0.6);
    DisplayObject->setHapticWorkspace(hduVector3Dd(-40, -40.0, -17.0), hduVector3Dd(95, 45.0, 17.0));
    // DisplayObject->enableLight(1);
    // DisplayObject->setLightAmbient(.6, .4, .4, 1.0);

    TriMesh* table = new TriMesh("Models/table1.obj");
    table->setRotation(hduVector3Dd(1.0, 0.0, 0.0), 0.0);
    table->setTranslation(hduVector3Dd(0.0, -50.0, 0.0));
    table->setTexture("Models/Textures/wood.jpg");
    table->setStiffness(0.9);
    table->setFriction(0.7, 0.5);
    table->setUnDraggable();
    DisplayObject->tell(table);

    DataTransportClass dataObject;  // Initialize an Object to transport data into the servoloop callback

    // Load cube1 model
    // dataObject.c1 = new TriMesh("Models/t21.obj");
    dataObject.c1 = new TriMesh("Models/box1.obj");
    // dataObject.c1 = new TriMesh("Models/pp11.obj");
    obj1 = dataObject.c1;
    dataObject.c1->setName("cube1");
    dataObject.c1->setShapeColor(1.0, 0.5, 0.65);
    dataObject.c1->setRotation(hduVector3Dd(1.0, 0.0, 0.0), 45.0);
    dataObject.c1->setTranslation(hduVector3Dd(0.0, 0.0, 0.0));
    // double x1 = (double)opti_obj_pos_init[0], y1 = (double)opti_obj_pos_init[1], z1 = (double)opti_obj_pos_init[2];
    // dataObject.c1->setTranslation(hduVector3Dd(x1, y1, z1));

    dataObject.c1->setStiffness(0.6);
    dataObject.c1->setDamping(0.1);
    dataObject.c1->setFriction(0.0, 0.9);
    dataObject.c1->setMass(0.9);
    DisplayObject->tell(dataObject.c1);  // Tell quickhaptics that cube exists

    // Load cube2 model
    // dataObject.c2 = new TriMesh("Models/c21.obj");
    dataObject.c2 = new TriMesh("Models/box1.obj");
    // dataObject.c2 = new TriMesh("Models/rr11.obj");
    obj2 = dataObject.c2;
    dataObject.c2->setName("cube2");
    dataObject.c2->setShapeColor(0.1, 0.5, 0.65);
    dataObject.c2->setRotation(hduVector3Dd(1.0, 0.0, 0.0), 45.0);
    dataObject.c2->setTranslation(hduVector3Dd(50.0, 0.0, 0.0));
    // double x2 = (double)hpti_obj_pos_init[0], y2 = (double)hpti_obj_pos_init[1], z2 = (double)hpti_obj_pos_init[2];
    // dataObject.c2->setTranslation(hduVector3Dd(x2, y2, z2));

    dataObject.c2->setStiffness(0.6);
    dataObject.c2->setDamping(0.1);
    dataObject.c2->setFriction(0.0, 0.9);
    dataObject.c2->setMass(0.9);
    DisplayObject->tell(dataObject.c2);  // Tell quickhaptics that cube exists

    // Load cube2-shadow model //////////////////////////////////////////////////////////////////////////
    // dataObject.c2_shadow = new TriMesh("Models/sp11.obj");
    dataObject.c2_shadow = new TriMesh("Models/box1.obj");
    obj2_shadow = dataObject.c2_shadow;
    dataObject.c2_shadow->setName("cube2_shadow");
    dataObject.c2_shadow->setShapeColor(1.0, 0.0, 0.0);
    dataObject.c2_shadow->setRotation(hduVector3Dd(1.0, 0.0, 0.0), 45.0);
    dataObject.c2_shadow->setTranslation(hduVector3Dd(50.0, 0.0, 0.0));
    dataObject.c2_shadow->setHapticVisibility(false);
    dataObject.c2_shadow->setGraphicVisibility(false);
    DisplayObject->tell(dataObject.c2_shadow);  // Tell quickhaptics that cube exists

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // dataObject.coll_sp = new TriMesh("Models/sp11.obj");
    dataObject.coll_sp = new TriMesh("Models/box1.obj");
    obj_coll_sp = dataObject.coll_sp;
    dataObject.coll_sp->setName("coll_sp");
    dataObject.coll_sp->setShapeColor(1.0, 0.0, 0.0);
    dataObject.coll_sp->setRotation(hduVector3Dd(1.0, 0.0, 0.0), 45.0);
    dataObject.coll_sp->setTranslation(hduVector3Dd(50.0, 0.0, 0.0));
    // dataObject.coll_sp->setScale(0.2);
    dataObject.coll_sp->setHapticVisibility(false);
    dataObject.coll_sp->setGraphicVisibility(false);
    DisplayObject->tell(dataObject.coll_sp);  // Tell quickhaptics that cube exists
    // //////////////////////////////////////////////////////////////////////////////////////////////////

    Text* text1 = new Text(18.0, " System Not Callibrated!", 0.15, 0.9);
    notCallibrate = text1;
    text1->setName("startupMsg");
    text1->setShapeColor(1.0, 0.05, 0.05);
    text1->setHapticVisibility(false);
    text1->setGraphicVisibility(true);
    DisplayObject->tell(text1);

    text1 = new Text(10.0, "Use right-click to push points (min 3 pts) then callibrate.", 0.1, 0.85);
    calibInstruction = text1;
    text1->setName("instructionMsg");
    text1->setShapeColor(1.0, 0.05, 0.05);
    text1->setHapticVisibility(false);
    text1->setGraphicVisibility(true);
    DisplayObject->tell(text1);

    Cursor* OmniCursor = new Cursor("Models/myCurser.obj");  // Load a cursor
    // Cursor* OmniCursor = new Cursor("Models/rr11.obj");  // Load a cursor
    TriMesh* cursorModel = OmniCursor->getTriMeshPointer();
    OmniCursor->setName("devCursor");  // Give it a name

    cursorModel->setShapeColor(0.35, 0.35, 0.35);
    // OmniCursor->scaleCursor(0.007);
    OmniCursor->scaleCursor(0.01);
    OmniCursor->setRelativeShapeOrientation(0.0, 0.0, 1.0, -90.0);
    // OmniCursor->setRelativeShapeOrientation(1.0, 0.0, 1.0, 90.0);

    //    OmniCursor->debugCursor(); //Use this function the view the location of the proxy inside the Cursor mesh
    DisplayObject->tell(OmniCursor);  // Tell QuickHaptics that the cursor exists

    DisplayObject->preDrawCallback(graphicsCallback);
    deviceSpace->startServoLoopCallback(startEffectCB, computeForceCB, stopEffectCB, &dataObject);  // Register the servoloop callback

    deviceSpace->button1UpCallback(button1UpCallback);

    // Create the GLUT menu
    glutCreateMenu(glutMenuFunction);
    // Add entries
    glutAddMenuEntry("push point", 0);
    glutAddMenuEntry("pop point", 1);
    glutAddMenuEntry("calibrate", 2);
    // glutAddMenuEntry("camera", 3);

    // attache to the menu
    glutAttachMenu(GLUT_RIGHT_BUTTON);

    // changing camera view
    DisplayObject->setDefaultCamera();

    float fov, nearplane, farplane;
    hduVector3Dd eyepoint, lookat, up;
    DisplayObject->getCamera(&fov, &nearplane, &farplane, &eyepoint, &lookat, &up);

    eyepoint[2] += 23.;  // pull back by 100
    // eyepoint[2] += eye[2];  // pull back by 100

    nearplane += 10.;  // recenter the haptic workspace (adjust by 20)
    farplane += 100.;
    DisplayObject->setCamera(fov + -15., nearplane, farplane, eyepoint, lookat, up);

    qhStart();  // Set everything in motion
}

void button1DownCallback(unsigned int ShapeID) {
    TriMesh* modelTouched = TriMesh::searchTriMesh(ShapeID);
    Box* buttonTouched = Box::searchBox(ShapeID);
}

void button1UpCallback(unsigned int ShapeID) {
    obj1->setHapticVisibility(true);
    obj2->setHapticVisibility(true);
}

void touchCallback(unsigned int ShapeID) {
    TriMesh* modelTouched = TriMesh::searchTriMesh(ShapeID);
}

void graphicsCallback() {
    /////////////////////////////////////////////////////////////////////////////////////////////// changing camera position
    //////////////////////////// getting screen tracker position
    Eigen::Vector3f p = (opti->rigidObjects)[screen_tracker]->position;

    // transform optitrack coordinate into haptic system
    Eigen::Vector3d sc_opti_to_hapt(0., 0., 0.);
    Eigen::Vector3d sc_opti_double;

    sc_opti_double[0] = (double)p[0];
    sc_opti_double[1] = (double)p[1];
    sc_opti_double[2] = (double)p[2];

    if (matrix_flag) {
#pragma opm parallel
        {
            t->transformPointFromSysAtoSysB(sc_opti_double, sc_opti_to_hapt);
            // std::cout << "----------Transformation Matrix----------" << std::endl;
            // std::cout << t->transformMat << std::endl;
            // std::cout << "----------Error----------" << std::endl;
            double ee = t->errorCalculation();
            std::cout << ">>>>>>>>>>>>>>>>>>>> Error: " << ee << std::endl;
        }
    }

    HDdouble screen_cam_x = (HDdouble)(sc_opti_to_hapt[0]);
    HDdouble screen_cam_y = (HDdouble)(sc_opti_to_hapt[1]);
    HDdouble screen_cam_z = (HDdouble)(sc_opti_to_hapt[2]);

    //////////////////////////// getting head tracker position
    HDdouble head_cam_x;
    HDdouble head_cam_y;
    HDdouble head_cam_z;
    if ((opti->rigidObjects)[head_tracker]) {
        Eigen::Vector3f p = (opti->rigidObjects)[head_tracker]->position;

        // transform optitrack coordinate into haptic system
        Eigen::Vector3d h_opti_to_hapt(0., 0., 0.);
        Eigen::Vector3d h_opti_double;

        h_opti_double[0] = (double)p[0];
        h_opti_double[1] = (double)p[1];
        h_opti_double[2] = (double)p[2];

        if (matrix_flag) {
#pragma opm parallel
            {
                t->transformPointFromSysAtoSysB(h_opti_double, h_opti_to_hapt);
                // std::cout << "----------Transformation Matrix----------" << std::endl;
                // std::cout << t->transformMat << std::endl;
                // std::cout << "----------Error----------" << std::endl;
                double ee = t->errorCalculation();
                std::cout << ">>>>>>>>>>>>>>>>>>>> Error: " << ee << std::endl;
            }
        }

        head_cam_x = (HDdouble)(h_opti_to_hapt[0]);
        head_cam_y = (HDdouble)(h_opti_to_hapt[1]);
        head_cam_z = (HDdouble)(h_opti_to_hapt[2]);

        screen_cam_x += head_cam_x;
        screen_cam_y += head_cam_y;
        screen_cam_z += head_cam_z;
    } else {
        if ((head_cam_x != NULL) && (head_cam_y != NULL) && (head_cam_z != NULL)) {
            screen_cam_x -= head_cam_x;
            screen_cam_y -= head_cam_y;
            screen_cam_z -= head_cam_z;
        }
        head_cam_x = 0.0;
        head_cam_y = 0.0;
        head_cam_z = 0.0;
    }

    //////////////////////////// updating camera position based on the screen tracker position
    QHGLUT* ds = QHGLUT::searchWindow("screen");
    // ds->setBackgroundColor(1.0, 0.0, 0.0);

    // changing camera view
    ds->setDefaultCamera();
    float fov, nearplane, farplane;
    hduVector3Dd eyepoint, lookat, up;
    ds->getCamera(&fov, &nearplane, &farplane, &eyepoint, &lookat, &up);

    eyepoint[0] += screen_cam_x;
    eyepoint[1] += screen_cam_y;
    eyepoint[2] += screen_cam_z;

    ds->setCamera(fov + -15., nearplane, farplane, eyepoint, lookat, up);

    /////////////////////////////////////////////////////////////////////////////////////////////// getting cursor position
    Cursor* localDeviceCursor = Cursor::searchCursor("devCursor");  // Get a pointer to the cursor
    hduVector3Dd localCursorPosition;
    localCursorPosition = localDeviceCursor->getPosition();  // Get the local cursor position in World Space
    // printf("--------------------------------------------------------- %f, %f, %f\n", localCursorPosition[0], localCursorPosition[1], localCursorPosition[2]);
    hcurr_temp[0] = (double)localCursorPosition[0];
    hcurr_temp[1] = (double)localCursorPosition[1];
    hcurr_temp[2] = (double)localCursorPosition[2];
    ///////////////////////////////////////////////////////////////////////////////////////////////
    if (matrix_flag) {
        notCallibrate->setGraphicVisibility(false);
        calibInstruction->setGraphicVisibility(false);
    }
    if (isOverlap) {
        obj2_shadow->setGraphicVisibility(true);
        obj2->setGraphicVisibility(false);
        if (isInside) {
            obj_coll_sp->setGraphicVisibility(true);
        }
    } else {
        obj2_shadow->setGraphicVisibility(false);
        obj2->setGraphicVisibility(true);
        obj_coll_sp->setGraphicVisibility(false);
    }
}

/***************************************************************************************
 Servo loop thread callback.  Computes a force effect. This callback defines the motion
 of the purple skull and calculates the force based on the "real-time" Proxy position
 in Device space.
****************************************************************************************/
void HLCALLBACK computeForceCB(HDdouble force[3], HLcache* cache, void* userdata) {
    DataTransportClass* localdataObject = (DataTransportClass*)userdata;  // Typecast the pointer passed in appropriately
    static int counter1 = 0;
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////// transforming object
    Eigen::Vector3f p = (opti->rigidObjects)[object_tracker]->position;

    ocurr_temp[0] = (double)p[0];
    ocurr_temp[1] = (double)p[1];
    ocurr_temp[2] = (double)p[2];

    // transform optitrack coordinate into haptic system
    Eigen::Vector3d opti_to_hapt(0., 0., 0.);
    Eigen::Vector3d opti_double;

    opti_double[0] = (double)p[0];
    opti_double[1] = (double)p[1];
    opti_double[2] = (double)p[2];

    if (matrix_flag) {
#pragma opm parallel
        {
            t->transformPointFromSysAtoSysB(opti_double, opti_to_hapt);
            // std::cout << "----------Transformation Matrix----------" << std::endl;
            // std::cout << t->transformMat << std::endl;
            // std::cout << "----------Error----------" << std::endl;
            double ee = t->errorCalculation();
            std::cout << ">>>>>>>>>>>>>>>>>>>> Error: " << ee << std::endl;
        }
    }

    HDdouble x = (HDdouble)(opti_to_hapt[0]);
    HDdouble y = (HDdouble)(opti_to_hapt[1]);
    HDdouble z = (HDdouble)(opti_to_hapt[2]);

    // printf("-------------------- %f, %f, %f\n", x, z, y);
    double rX;
    double rY;
    double rZ;

    hduMatrix rZM;
    hduMatrix rYM;
    hduMatrix rXM;

    hduMatrix rxyz;

    // #pragma opm parallel
    //     {
    rX = (opti->rigidObjects)[1]->rotation[0] * 3.14159 / 180.0;
    rY = (opti->rigidObjects)[1]->rotation[1] * 3.14159 / 180.0;
    rZ = (opti->rigidObjects)[1]->rotation[2] * 3.14159 / 180.0;

    rZM = hduMatrix::createRotationAroundZ(rZ);
    rYM = hduMatrix::createRotationAroundY(rY);
    rXM = hduMatrix::createRotationAroundX(rX);
    rxyz = rZM * rYM * rXM;
    // }

    localdataObject->c1->setTransform(rxyz);       // rotate the skull with the optitrack trackerocurr_temp
    localdataObject->c1->setTranslation(x, y, z);  // move the skull with the optitrack tracker

    // localdataObject->c1->setScaleInPlace(0.3);

    double rbX;
    double rbY;
    double rbZ;
#pragma opm parallel
    {
        //////////////////////////////////////////////////////////// opti to Quatrenion (nut)
        Eigen::Vector3d rotation(rX, rY, rZ);
        double angle = rotation.norm();
        Eigen::Vector3d axis = rotation.normalized();
        Eigen::Quaterniond q(Eigen::AngleAxisd(angle, axis));

        opti_obj_q.x = q.x();
        opti_obj_q.y = q.y();
        opti_obj_q.z = q.z();
        opti_obj_q.w = q.w();

        // reactphysics3d::Vector3 new_opti_obj_pos(opti_to_hapt[0], opti_to_hapt[1], opti_to_hapt[2]);
        // updated_opti_obj_pos = new_opti_obj_pos;
        updated_opti_obj_pos[0] = opti_to_hapt[0];
        updated_opti_obj_pos[1] = opti_to_hapt[1];
        updated_opti_obj_pos[2] = opti_to_hapt[2];
        //////////////////////////////////////////////////////////// hapti to Quatrenion (bolt)
        /*
        getting rotation angles from rotation matrix:

        R = |r00, r01, r02|
            |r10, r11, r12|
            |r20, r21, r22|

        rx = atan2(r21, r22)
        ry = atan2(-r20, sqrt(pow(r21, 2) + pow(r22, 2)))
        rz = atan2(r10, r00)
        */

        hduMatrix rmat = localdataObject->c2->getRotation();

        rbX = std::atan2(rmat[2][1], rmat[2][2]) * 180.0 / 3.14159;
        rbY = std::atan2(-rmat[2][0], std::sqrt(rmat[2][1] * rmat[2][1] + rmat[2][2] * rmat[2][2])) * 180.0 / 3.14159;
        rbZ = std::atan2(rmat[1][0], rmat[0][0]) * 180.0 / 3.14159;

        Eigen::Vector3d rotationB(rbX, rbY, rbZ);
        double angleBolt = rotationB.norm();
        Eigen::Vector3d axisBolt = rotationB.normalized();
        Eigen::Quaterniond qB(Eigen::AngleAxisd(angleBolt, axisBolt));

        hpti_obj_q.x = qB.x();
        hpti_obj_q.y = qB.y();
        hpti_obj_q.z = qB.z();
        hpti_obj_q.w = qB.w();

        hduVector3Dd h_pos = localdataObject->c2->getTranslation();
        double hx = h_pos[0], hy = h_pos[1], hz = h_pos[2];

        updated_hpti_obj_pos[0] = hx;
        updated_hpti_obj_pos[1] = hy;
        updated_hpti_obj_pos[2] = hz;

        // setting the translation of the shadow object to the actual haptic object
        localdataObject->c2_shadow->setTransform(rmat);
        localdataObject->c2_shadow->setTranslation(h_pos);

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // #pragma opm parallel
        //         {
        //             if (isOverlap) {
        //                 for (int i = 0; i < 482; i++) {
        //                     hduVector3Dd v = localdataObject->c2->getVertex(i);
        //                     reactphysics3d::Vector3 temp_pos;
        //                     temp_pos.x = v[0];
        //                     temp_pos.y = v[1];
        //                     temp_pos.z = v[2];
        //                     isInside = body1_check->testPointInside(temp_pos);
        //                     if (isInside) {
        //                         // localdataObject->coll_sp->setTransform(rmat);
        //                         localdataObject->coll_sp->setTranslation(v);
        //                         std::cout << "-------------------------------------------inside!" << std::endl;
        //                     }
        //                 }
        //             }
        //         }

        double dist = std::sqrt(std::pow((opti_to_hapt[0] - hx), 2) + std::pow((opti_to_hapt[1] - hy), 2) + std::pow((opti_to_hapt[2] - hz), 2));
        // std::cout << "-----------------------dist: " << dist << std::endl;
    }
    // std::cout << "====================================O" << updated_opti_obj_pos[0] << " " << updated_opti_obj_pos[1] << " " << updated_opti_obj_pos[2] << std::endl;
    // std::cout << "------------------------------------H" << updated_hpti_obj_pos[0] << " " << updated_hpti_obj_pos[1] << " " << updated_hpti_obj_pos[2] << std::endl;
    // std::cout << "====================================OQ" << opti_obj_q.x << " " << opti_obj_q.y << " " << opti_obj_q.z << " " << opti_obj_q.w << std::endl;
    // std::cout << "====================================HQ" << hpti_obj_q.x << " " << hpti_obj_q.y << " " << hpti_obj_q.z << " " << hpti_obj_q.w << std::endl;

    if (isOverlap)
        std::cout << "=================================Overlap!" << std::endl;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#pragma opm parallel
    {
        HDdouble nominalMaxContinuousForce;
        hdGetDoublev(HD_NOMINAL_MAX_CONTINUOUS_FORCE, &nominalMaxContinuousForce);
        if (isOverlap) {
            f_decay_count += 0.000005;
            force[0] = nominalMaxContinuousForce * 1.0 / f_decay_count;
            force[1] = nominalMaxContinuousForce * 1.0 / f_decay_count;
            force[2] = nominalMaxContinuousForce * 1.0 / f_decay_count;
            if (f_decay_count > 3000) f_decay_count = 1.0;
        } else {
            force[0] = 0.0;
            force[1] = 0.0;
            force[2] = 0.0;
            f_decay_count = 1.0;
        }
    }
}

/******************************************************************************
 Servo loop thread callback called when the effect is started.
******************************************************************************/
void HLCALLBACK startEffectCB(HLcache* cache, void* userdata) {
    DataTransportClass* localdataObject = (DataTransportClass*)userdata;
    printf("Custom effect started\n");
}

/******************************************************************************
 Servo loop thread callback called when the effect is stopped.
******************************************************************************/
void HLCALLBACK stopEffectCB(HLcache* cache, void* userdata) {
    printf("Custom effect stopped\n");
}

/*******************************************************************************
 Given the position of the two charges in space,
 calculates the (modified) coulomb force.
*******************************************************************************/
hduVector3Dd forceField(hduVector3Dd Pos1, hduVector3Dd Pos2, HDdouble Multiplier, HLdouble Radius) {
    hduVector3Dd diffVec = Pos2 - Pos1;  // Find the difference in position
    double dist = 0.0;
    hduVector3Dd forceVec(0, 0, 0);

    HDdouble nominalMaxContinuousForce;
    hdGetDoublev(HD_NOMINAL_MAX_CONTINUOUS_FORCE, &nominalMaxContinuousForce);  // Find the max continuous for that the device is capable of

    dist = diffVec.magnitude();

    if (dist < Radius * 2.0)  // Spring force (when the model and cursor are within a 'sphere of influence'
    {
        diffVec.normalize();
        forceVec = (Multiplier)*diffVec * dist / (4.0 * Radius * Radius);
        static int i = 0;
    } else  // Inverse square attraction
    {
        forceVec = Multiplier * diffVec / (dist * dist);
    }

    for (int i = 0; i < 3; i++)  // Limit force calculated to Max continuouis. This is very important because force values exceeding this value can damage the device motors.
    {
        if (forceVec[i] > nominalMaxContinuousForce)
            forceVec[i] = nominalMaxContinuousForce;

        if (forceVec[i] < -nominalMaxContinuousForce)
            forceVec[i] = -nominalMaxContinuousForce;
    }

    return forceVec;
}

void glutMenuFunction(int MenuID) {
    if (MenuID == 0) {
        haptic_ptsB.push_back(hcurr_temp);
        optitrack_ptsA.push_back(ocurr_temp);
    }

    if (MenuID == 1) {
        if (haptic_ptsB.size() >= 0)
            haptic_ptsB.pop_back();
        if (optitrack_ptsA.size() >= 0)
            optitrack_ptsA.pop_back();
    }

    if (MenuID == 2) {
        // calibrate
        t = std::make_shared<coordinateTransform>(optitrack_ptsA, haptic_ptsB);
        t->calculateTransformMatrix();
        matrix_flag = true;
    }

    // if (MenuID == 3) {
    //     eye[0] += 3.;
    //     eye[1] += 1.;
    //     eye[2] += 5.;
    // }
}
