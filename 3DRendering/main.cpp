#include "raylib.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <sstream>
#include <thread>
#include <cmath>
#include <array>
#include <chrono>
#include "rlgl.h"
#include "../consts.hpp"
using namespace std;

const long TIME_STEP = pow(10,9)*loopPerSample*dt; // 100000000; //ns
const float SPEED_UP = 1;

Vector3 vectToVect3(vector<float> vect) {
    //cout << "Proj [x,z,y]" + to_string(vect[0]) + ", " + to_string(vect[2]) + ", " + to_string(vect[1]) << "\n"; 
    return {vect[0],vect[2],vect[1]}; 
 }
float dotProduct(vector<float> vect_A, vector<float> vect_B)
{
 
    float product = 0;
 
    for (int i = 0; i < (int)vect_A.size(); i++) {
        product = product + vect_A[i] * vect_B[i];
    }
    return product;
}

void normalize(std::vector<float>* arr, int scale) {
    
    double mod = 0.0;

    for (size_t i = 0; i < (*arr).size(); ++i) {
        mod += (*arr)[i] * (*arr)[i];
    }

    double mag = sqrt(mod);

    if (mag == 0) {
        throw std::logic_error("The input vector is a zero vector");
    }

    for (size_t i = 0; i < (*arr).size(); ++i) {
        (*arr)[i] = ((*arr)[i] / mag)*scale;
    }
}
float magnitude(std::vector<float>* arr) {
    float sum = 0; 
    for (float val: *arr) {
        sum += pow(val,2);
    }
    sum = sqrt(sum);
    return sum;
}
vector<float> projVect(vector<float> vect, vector<float> normal) {
    vector<float> proj;
    vector<float> reject;
    normalize(&normal, 1);
    float dot = dotProduct(vect, normal);
    for (int i = 0; i < 3; i++) {
        proj.push_back(vect[i] - normal[i]*dot);
    }
    //cout << "\nVect Projection";
    /*cout << "\n";
    // cout << "\nVector";
    for (int i = 0; i < 3; i++) {
        cout << vect[i] << " ";
    }
    cout << "\n";
    */
    return proj; 
}
void log_vec(std::vector<float>* arr, int scale) {
    for (size_t i = 0; i < (*arr).size(); ++i) {
        (*arr)[i] = log((*arr)[i]);
    }
    
}

int main(void)
{
    ifstream programCode ("../results.txt");
    vector<vector<float>> pos;
    vector<vector<float>> ang;
    vector<vector<float>> vel;
    vector<vector<float>> acc;
    vector<float> timeElapsed;
    vector<float> fuel;
    vector<float> normal = {0.0f, 0.0f, 1.0f};
    int scale = 10;
    string data;
    
    if (!programCode.is_open())
    {
        cout << "Can't open file" << "\n";
        exit(1);
    }
    while (getline(programCode, data)) {
        stringstream ss(data);
        // TIME,POS[X],POS[Y],POS[Z],Ang[X],Ang[Y],Ang[Z],fuel
        
        if (ss.good()) {
            string substr;
            getline(ss, substr, ',');
            //cout << stof(substr) << "\n";
            timeElapsed.push_back(stof(substr));

            vector<float> pos_arr;
            getline(ss, substr, ',');
            pos_arr.push_back(stof(substr)); 
            
            getline(ss, substr, ',');
            pos_arr.push_back(stof(substr));

            getline(ss, substr, ',');
            pos_arr.push_back(stof(substr));

            pos.push_back(pos_arr);
     
            vector<float> ang_arr;
            getline(ss, substr, ',');
            ang_arr.push_back(stof(substr));

            getline(ss, substr, ',');
            ang_arr.push_back(stof(substr));

            getline(ss, substr, ',');
            ang_arr.push_back(stof(substr));

            ang.push_back(ang_arr);

            vector<float> vel_arr;
            getline(ss, substr, ',');
            vel_arr.push_back(stof(substr));

            getline(ss, substr, ',');
            vel_arr.push_back(stof(substr));

            getline(ss, substr, ',');
            vel_arr.push_back(stof(substr));

            vel.push_back(vel_arr);

            vector<float> acc_arr;
            getline(ss, substr, ',');
            acc_arr.push_back(stof(substr));

            getline(ss, substr, ',');
            acc_arr.push_back(stof(substr));

            getline(ss, substr, ',');
            acc_arr.push_back(stof(substr));

            acc.push_back(acc_arr);

            

            getline(ss, substr);
            fuel.push_back(stof(substr));
            
        }
        
    }
    
    // Initialization
    //--------------------------------------------------------------------------------------
    const int screenWidth = 2000;
    const int screenHeight = 1000;

    InitWindow(screenWidth, screenHeight, "Rocket Simulation Render");

    // Define the camera to look into our 3d world
    Camera3D camera = { 0 };
    camera.position = (Vector3){ 0.0f, 17.0f, 17.0f };  // Camera position
    camera.target = (Vector3){ 0.0f, 0.0f, 0.0f };      // Camera looking at point
    camera.up = (Vector3){ 0.0f, 1.0f, 0.0f };          // Camera up vector (rotation towards target)
    camera.fovy = 45.0f;                                // Camera field-of-view Y
    camera.projection = CAMERA_PERSPECTIVE;             // Camera mode type

    Model rocket = LoadModel("./Assets/rocketNeworigin.obj"); 
    
    Vector3 cubePosition = { 0.0f, 0.0f, 0.0f };
    //DrawModel(rocket, cubePosition, 1.0f, RED);
    //auto timeNow = std::chrono::high_resolution_clock::now();
    //auto timePrevious = timeNow;
    
    //float simTime = timeElapsed.front();
    //float simTimePast = 0.0f;

    //SetTargetFPS(10);               // Set our game to run at 60 frames-per-second
    //--------------------------------------------------------------------------------------
    string timeELapsed_str = "Time Elapsed (s) " + std::to_string(timeElapsed.front());
    //simTime = timeElapsed.front();
    const char* timeElapsed_charPtr = timeELapsed_str.c_str();
    timeElapsed.erase(timeElapsed.begin());

    float ang_x = ang.front()[0];
    float ang_y = ang.front()[1];
    float ang_z = ang.front()[2];
    ang.erase(ang.begin());

    vector<float> velVec = vel.front();
    float vel_mag = magnitude(&velVec);
    string vel_str = "Vel: \nX:" + to_string(velVec[0]) + "\nY:" + to_string(velVec[1])  + "\nZ:" + to_string(velVec[2]) + "\nMag:" + to_string(vel_mag); 
    const char* vel_ptr = vel_str.c_str();
    normalize(&(velVec),scale);
    vector<float> vel_proj = projVect(velVec, normal);
    float vel_x = velVec[0];
    float vel_y = velVec[1];
    float vel_z = velVec[2];
    vel.erase(vel.begin());

    vector<float> accVec = acc.front();
    float acc_mag = magnitude(&accVec);
    string acc_str = "Acc: \nX:" + to_string(accVec[0]) + "\nY:" + to_string(accVec[1])  + "\nZ:" + to_string(accVec[2]) + "\nMag:" + to_string(acc_mag);; 
    const char* acc_ptr = acc_str.c_str();
    normalize(&(accVec),scale);
    vector<float> acc_proj = projVect(accVec, normal);
    float acc_x = accVec[0];
    float acc_y = accVec[1];
    float acc_z = accVec[2];
    acc.erase(acc.begin());

    string fuel_str = "Fuel (kg) " + std::to_string(fuel.front());
    const char* fuel_ptr = fuel_str.c_str();
    fuel.erase(fuel.begin());

    auto timePrevious = std::chrono::high_resolution_clock::now();
    auto timeNow = std::chrono::high_resolution_clock::now();
    long long dlta = 0;
    // Main game loop
    while (!WindowShouldClose())    // Detect window close button or ESC key
    {
        // Update
        //----------------------------------------------------------------------------------
        // TODO: Update your variables here
        //----------------------------------------------------------------------------------
        /*
        float wait_time = (simTime - simTimePast) - pow(10,-9)*std::chrono::duration_cast<std::chrono::nanoseconds>(timeNow-timePrevious).count();
        cout << (wait_time*1000ms).count() << "\n";
        cout << wait_time << "\n";
        simTimePast = simTime;
        timePrevious = std::chrono::high_resolution_clock::now();
        

        std::this_thread::sleep_for((wait_time*1000ms));
        timeNow = std::chrono::high_resolution_clock::now();
        */
        timePrevious = timeNow;
        timeNow = std::chrono::high_resolution_clock::now(); 
        //cout << std::chrono::duration<float>(timeNow-timePrevious).count() << '\n';
        dlta += std::chrono::duration_cast<std::chrono::nanoseconds>(timeNow-timePrevious).count();
        if (dlta > TIME_STEP/SPEED_UP) {
            dlta -= TIME_STEP/SPEED_UP;
            timePrevious = timeNow;
            timeELapsed_str = "Time Elapsed (s) " + std::to_string(timeElapsed.front());
            //simTime = timeElapsed.front();
            timeElapsed_charPtr = timeELapsed_str.c_str();
            timeElapsed.erase(timeElapsed.begin());

            ang_x = ang.front()[0];
            ang_y = ang.front()[1];
            ang_z = ang.front()[2];
            ang.erase(ang.begin());

            velVec = vel.front();
            vel_mag = magnitude(&velVec);
            vel_str = "Vel: \nX:" + to_string(velVec[0]) + "\nY:" + to_string(velVec[1])  + "\nZ:" + to_string(velVec[2]) + "\nMag:" + to_string(vel_mag); 
            vel_ptr = vel_str.c_str();
            normalize(&(velVec),scale);
            vel_proj = projVect(velVec, normal);
            vel_x = velVec[0];
            vel_y = velVec[1];
            vel_z = velVec[2];
            vel.erase(vel.begin());

            accVec = acc.front();
            acc_mag = magnitude(&accVec);
            acc_str = "Acc: \nX:" + to_string(accVec[0]) + "\nY:" + to_string(accVec[1])  + "\nZ:" + to_string(accVec[2]) + "\nMag:" + to_string(acc_mag);; 
            acc_ptr = acc_str.c_str();
            normalize(&(accVec),scale);
            acc_proj = projVect(accVec, normal);
            acc_x = accVec[0];
            acc_y = accVec[1];
            acc_z = accVec[2];
            acc.erase(acc.begin());

            fuel_str = "Fuel (kg) " + std::to_string(fuel.front());
            fuel_ptr = fuel_str.c_str();
            fuel.erase(fuel.begin());
        }
        // Draw
        //----------------------------------------------------------------------------------
        BeginDrawing();

            ClearBackground(RAYWHITE);

            BeginMode3D(camera);
                
                rlPushMatrix();
                // (x,z,y) counter clockwise is positive angle (degrees) 

                //rlTranslatef(0,) //need to center on pivot point (bottom center of rocket)
                rlRotatef(ang_x,1,0,0);
                rlRotatef(ang_y,0,0,1);
                rlRotatef(ang_z,0,1,0);
                //rlRotatef(90,0,1,0); //

                DrawModel(rocket, cubePosition, .005f, WHITE);
                rlPopMatrix();
                
                DrawLine3D(cubePosition, {vel_x,vel_z,vel_y}, RED);
                //cout << "Vel: " + to_string(vel_x) + "  " + to_string(vel_z) + " " + to_string(vel_y) << "\n";
                DrawLine3D(cubePosition, vectToVect3(vel_proj), BLACK);
                DrawLine3D(cubePosition, {acc_x,acc_z,acc_y}, BLUE);
                DrawLine3D(cubePosition, vectToVect3(acc_proj), BLACK);
                DrawCubeWires(cubePosition, 2.0f, 2.0f, 2.0f, MAROON);
                DrawGrid(20, 1.0f);

            EndMode3D();
            //string timeELapsed_str = "Time Elapsed (s) " + std::to_string(timeElapsed.front());
            //const char* timeElapsed_charPtr = timeELapsed_str.c_str();
            //timeElapsed.erase(timeElapsed.begin());

            DrawText(timeElapsed_charPtr, 10, 40, 20, DARKGRAY);
            DrawText(fuel_ptr, 60, 80, 20, DARKGRAY);
            DrawText(vel_ptr, 60, 120, 20, DARKGRAY);
            DrawText(acc_ptr, 60, 280, 20, DARKGRAY);
            DrawText("-X", 440, 490, 20, DARKGRAY);
            DrawText("-Y", 950, 190, 20, DARKGRAY);

            DrawFPS(10, 10);

        EndDrawing();

        //timeNow = std::chrono::high_resolution_clock::now();
        //----------------------------------------------------------------------------------
    }
    UnloadModel(rocket);
    // De-Initialization
    //--------------------------------------------------------------------------------------
    CloseWindow();        // Close window and OpenGL context
    //--------------------------------------------------------------------------------------
    
    return 0;
}
