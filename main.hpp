#include "consts.hpp"
const double PI = 3.14159;
using namespace std;
class Tank;
class Rocket;
class Dynamics;
class Environment;
class Control;
class Shape {
  public:
    Shape(float SA[], double COV[], float origin[], float volume) 
    : sa(SA), cov(COV), P(origin), vol(volume) {
      updateCOVLocal();
    }
    float* getSA();
    double* getI();
    float* getOrigin();
    double* updateCOVLocal();
    double* getCOVLocal();
    double* getCOV();
    virtual double* updateI(float);
    virtual void updateCOM(float*, float, float, float);
    virtual void updateFilledVolume(Tank* t, double dltaM);
    virtual double* updateIFuel(Tank* t);
    
  protected:
    float* const sa; //surface area
    double* const cov; // center of volume (dry mass not considering liquid that is handled by Tank class)
    float* const P; // distance from P to local P
    double* const covLocal = new double[DIM];
    const float vol;
    double* I = new double[DIM];
};
class RectagularPrism : public Shape
{
  public:
    RectagularPrism(float SA[], double COV[], float P[], double volume, float leng, float wid, float hei) 
    : Shape(SA, COV, P, volume), length(leng), width(wid), height(hei) {}
    float getAreaLW();
    float getLength();
    float getWidth();
    float getHeight();
    void updateFilledVolume(Tank* t, double dltaM) override;
    double* updateIFuel(Tank* t) override;
    double* updateI(float mass) override;
    void updateCOM(float* com, float fuelH, float fuelM, float dry_mass) override;
    
  private:
    float length;
    float width;
    float height;
};
class Sphere : public Shape
{
   public:
    Sphere(float SA[], double COV[], float dist[], double volume, float rad) 
    : Shape(SA, COV, dist, volume), radius(rad) {}
    float getRadius();
    void updateFilledVolume(Tank* t, double dltaM) override;
    float getVolume();
    double* updateIFuel(Tank* t) override;
    //hollow sphere
    double* updateI(float mass) override;
    void updateCOM(float* com, float fuelH, float fuelM, float dry_mass) override;
  private:
    float radius;
};

/*
 * currently we are assuming that the fuel stays perpendicular to the z axis of the rocket which obviously will need to changed in the future.
 */
class Tank {
  public:
    Tank(Shape* s, float dry_mass, float massF, double mass_vel_out, double fuelD, float volFill, float fuelHeight) 
    : s(s), dry_mass(dry_mass), fuelM(massF), mass_vel_out(mass_vel_out), volumeFilled(volFill), densityFuel(fuelD), fuelH(fuelHeight) 
    {
      s->updateI(dry_mass);
      s->updateCOM(com, fuelH, fuelM, dry_mass);
    }
    //calculate momentum of inertia caused by fuel, calculate base on container of fuel
    double* updateIFuel(Shape* s);
    double* updateIFuelHelper(RectagularPrism* s);
    //http://blitiri.blogspot.com/2014/05/mass-moment-of-inertia-of-hemisphere.html
    double* updateIFuelHelper(Sphere* s);
    //change in mass is in kg 
    int updateMass(double dt);
    float getDryMass();
    float getFuelMass();
    float getFuelHeight();
    float getMass();
    double* getIFuel();
    float* getCOM();
    Shape* getShape();
    void updateFilledVolume(Shape* s, double dltaM);
    void updateFilledVolumeHelper(RectagularPrism* s, double dltaM);
    void updateFilledVolumeHelper(Sphere* s, double dltaM);

    
    

    
    

  private:
    Shape* s; //simple shape with known area, volume, other parameters, etc
    const float dry_mass;
    float fuelH;
    const float densityFuel;
    float volumeFilled;
    float fuelM;
    double* const fuelI = new double[DIM];
    float* const com = new float[DIM]; // center of mass (meters) relative to P. this only accounts for wet mass
    double mass_vel_out; //kg/s
};  

class Dynamics { 
public:
  class State {
      public:
        State(double* angle, double* position, double* velo, double* ang_velo, double* acce, double* ang_acce, int m)
        {
          ang = angle; //euler angle (this is not the same as inertial frame)
          pos = position; //distance from global origin to rocket origin
          vel = velo;
          ang_vel = ang_velo; //rotation in inertial frame
          ang_acc = ang_acce; //rotation in inertial frame
          acc = acce;
          mass = m;
          timeElapsed = 0;
          prev = new State(*this);
        }
        double* ang; //global frame variables
        double* pos;
        double* vel;
        double* ang_vel;
        double* acc;
        double* ang_acc;
        double timeElapsed;
        int mass;
        State* prev;

        void updatePrev(State* curr);
        string to_string(); 
  };
  int const NUM_FORCES = 1; //change this value and array of forces and moments
  Dynamics (Rocket* rocket, Environment* env, Control* control, State* state) {
    r = rocket;
    e = env;
    c = control;
    s = state;
     //inertial force creating a moment

    //List of forces on the rocket
    fDrag = new double[DIM]; // just update forces so we don't have to make new memory have time the function is called
    fDrag_Ang = new double[DIM]; //drag caused by angular motion
    fThrust = new double[DIM];
    fG = new double[DIM];
    inerM = new double[DIM];

    //List of moment on the rocket;
    mDrag = new double[DIM]; // just update forces so we don't have to make new memory have time the function is called
    mThrust = new double[DIM];
    mG = new double[DIM];


    netF = new double[DIM]; //net force
    netM = new double[DIM]; //net moment
    
  }
  string to_string();
  string pos_string(double elaspedTime);
  void RK4(double dt);
  void updateState(double dt);
  void euler(double dt);

  double dPosdt(double dt, double vel, double acc);
  void crossProduct(double dist[], double force[], double moment[]);
  /*
   * @param double[] air, model air resistance as a vector 
   * vel is the velocity of the rocket, need to find relative wind vector
   * cd is coffiecient of drag of the rocket in each vector direction / it is proportional to surface area of that part of the rocket
   * TODO account for air resistance that flow around the rocket
   */
  void updateFdrag();
  /*
   * the moment created from the force of drag
   */
  void updateMdrag();



  void updateFg();
  void updateMg();

   /*
   * given by the control block
   */
  void updateFThrust(Dynamics::State* s, double dt);
  void updateMThrust();
  //moment caused by movement
  void updateInerMom();

  void updateNetF(double** forces);
  void updateNetM(double** moments);

  void globalToRocketFrame(double* loc, double* glo);
  void rocketToGlobalFrame(double* loc, double* glo);
  
  void clearNet();
  void updateAcc();
  void updateVel(double dt);
  double* getVel();
  double* getAcc();
  double* getAngVel();
  double* getAngAcc();
  void updatePos(double dt);

  void updateVel(double dt, double acc[]);
  void updatePos(double dt, double vel[]);

  void updateAngAcc();
  void updateAngVel(double dt);
  void updateAng(double dt);
  void updateAng(double dt, double ang_vel[]);

  //vector<boost::tuple<double, double>> time_vs_height;
  State* s;
private:
  
  
  double* inerM; //inertial force creating a moment


  //List of forces on the rocket
  double* fDrag; // just update forces so we don't have to make new memory have time the function is called
  double* fDrag_Ang;
  double* fThrust;
  double* fG;

  //List of moment on the rocket;
  double* mDrag; // just update forces so we don't have to make new memory have time the function is called
  double* mThrust;
  double* mG;


  double* netF; //net force
  double* netM; //net moment
  Rocket* r;
  Environment* e;
  Control* c;
  
};
class Control {
  public:
    Control(double* t) {
      thrust = t;
      initProfile();
    }
    double* updateThrust(Dynamics::State* s, Rocket* r, float fuelMass, double dt);
    double* getThrust();
    double** initProfile();
    void updateProfile(double dt);
    private:
      double* thrust;
      double** flight_profile;
      int time_itr = 0; //increments based on if enough time elapsed to go to next step
      double secs = 0; //time elapsed, resets after each iteration of time_itr
};
class Rocket {
  public:
    Rocket(Shape** shapes, double* centerofVol, float coeffDrag [], int d_m, float surfaceArea [], float m_vel_out, Tank** t, int numTank, int numShape, int* mParts, float max_t) 
      : s(shapes), cov(centerofVol), cd(coeffDrag), dry_mass(d_m), sa(surfaceArea), mass_vel_out(m_vel_out), massParts(mParts), tanks(t), 
        numTanks(numTank), numShapes(numShape), max_thrust(max_t)
        { 
          updateCOP();
          updateMass(0);
          updateIDry();
          updateI();
          updateComDry();
          updateCOM();
        }

    
    /*
    * updatePARAM vs getPARAM (sometimes you just want to get the PARAM vs getting and updating it) this is done to avoid recalculating values that we already know
    *
    */
    void updateComDry();
    Tank** getTanks();
    double* updateIDry();
    double* updateI();
    double* updateCOM();
    double* updateCOP();
    double* getI();
    int updateMass(double dt);
    /*
    * return Coeffient of drag based on the geometry of the surface of each direction of the rocket
    */
    float* getCD();
    double* getCOM();
    double* getCOP();
    float* getSA();
    int getMass();
    const int* getP();
    const double getMaxThrust();
    const int getNumParts();
    const int getNumTanks();
  private:
    //constants 

    Shape** const s; //array of shapes that describe the rocket
    const int* massParts; //parallel array with s, mass of each part
    const int numShapes;

    double* const cov; // center of volume, measure from P
    float* const cd;  //coefficient of drag
    float* const sa; //surface area
    double* const I = new double[DIM]; //moment of Inertia
    double* const IDry = new double[DIM];
    const int P[DIM] = {0,0,0}; //gimbal base, all rotations are taken with respect to this point, distance from rocket origin 
    const int dry_mass; 
    
    const float mass_vel_out; 
    const float max_thrust;
    // variables



    int tot_mass; // Newtons , tot_mass = dry_mass + m_tank1 + m_tank2 + ... + m_tankN
    double* const com = new double[DIM]; //center of mass
    double* const comDry = new double[DIM]; //center of dry mass constant
    double* const cop = new double[DIM]; //center of pressure
    
    Tank** tanks; //arbitary number of fuel tanks (lander should only have two but future designs may have different numbers)
    const int numTanks;
};