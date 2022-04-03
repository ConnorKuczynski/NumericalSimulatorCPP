#include "raylib.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <sstream>
#include "rlgl.h"
using namespace std;
int main(void)
{
    ifstream programCode ("../results.txt");
    vector<vector<float>> pos;
    vector<vector<float>> ang;
    vector<float> timeElapsed;
    vector<float> fuel;

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
            //cout << pos.front()[0] << "\n";
            vector<float> ang_arr;
            getline(ss, substr, ',');
            ang_arr.push_back(stof(substr));

            getline(ss, substr, ',');
            ang_arr.push_back(stof(substr));

            getline(ss, substr, ',');
            ang_arr.push_back(stof(substr));

            ang.push_back(ang_arr);

            getline(ss, substr);
            fuel.push_back(stof(substr));
            
        }
        
    }
    
    // Initialization
    //--------------------------------------------------------------------------------------
    const int screenWidth = 2000;
    const int screenHeight = 1000;

    InitWindow(screenWidth, screenHeight, "raylib [core] example - 3d camera mode");

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
    

    SetTargetFPS(60);               // Set our game to run at 60 frames-per-second
    //--------------------------------------------------------------------------------------

    // Main game loop
    while (!WindowShouldClose())    // Detect window close button or ESC key
    {
        // Update
        //----------------------------------------------------------------------------------
        // TODO: Update your variables here
        //----------------------------------------------------------------------------------
        
        string timeELapsed_str = "Time Elapsed (s) " + std::to_string(timeElapsed.front());
        const char* timeElapsed_charPtr = timeELapsed_str.c_str();
        timeElapsed.erase(timeElapsed.begin());

        int ang_x = ang.front()[0];
        int ang_y = ang.front()[1];
        int ang_z = ang.front()[2];
        ang.erase(ang.begin());

        string fuel_str = "Fuel (kg) " + std::to_string(fuel.front());
        const char* fuel_ptr = fuel_str.c_str();
        fuel.erase(fuel.begin());
        // Draw
        //----------------------------------------------------------------------------------
        BeginDrawing();

            ClearBackground(RAYWHITE);

            BeginMode3D(camera);
                
                rlPushMatrix();
                //rlTranslatef(0,) //need to center on pivot point (bottom center of rocket)
                rlRotatef(ang_x,1,0,0);
                rlRotatef(ang_y,0,1,0);
                rlRotatef(ang_z,0,0,1);
                DrawModel(rocket, cubePosition, .005f, BLACK);
                rlPopMatrix();
                DrawCubeWires(cubePosition, 2.0f, 2.0f, 2.0f, MAROON);

                DrawGrid(10, 1.0f);

            EndMode3D();
            //string timeELapsed_str = "Time Elapsed (s) " + std::to_string(timeElapsed.front());
            //const char* timeElapsed_charPtr = timeELapsed_str.c_str();
            //timeElapsed.erase(timeElapsed.begin());

            DrawText(timeElapsed_charPtr, 10, 40, 20, DARKGRAY);
            DrawText(fuel_ptr, 60, 80, 20, DARKGRAY);

            DrawFPS(10, 10);

        EndDrawing();
        //----------------------------------------------------------------------------------
    }
    UnloadModel(rocket);
    // De-Initialization
    //--------------------------------------------------------------------------------------
    CloseWindow();        // Close window and OpenGL context
    //--------------------------------------------------------------------------------------
    
    return 0;
}
