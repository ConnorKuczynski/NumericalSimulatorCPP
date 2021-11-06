#include <iostream>
#include <Eigen/Dense>

const int DIM = 2; // simulator dimension
const int X = 0; //enumerate accessing values from pos/vel/acc arrays 
const int Y = 1;
const int Z = 2;

/*
 * stores all the physical parameters of the rocket
 * list of param:
 * dry mass, fuel mass, oxidizer mass, length, width, height, (idea is to model each rocket part as a simple shape where we know the moment of inertia
 *       For example, we seperate the fuel tanks from the dry body, and model them as a sphere (or the best approximate shape), while body is a rectangle)
 * 
 * Summing the center of mass of each individual shape divided by total mass is equal to the center of mass of entire rocket 
 * 
 * we will model the initial rocket as a rectangular prism with a moment of inertia equal to 
 * 
 */
class Tank {
public:
  Tank() {
  }
  int updateMass(float dt)
  {
    dltaM = mass_vel_out*dt;

    updateCOM(calcDltaH(dltaM));

  }
  /*
   *  we only need to calc z because h_fuel changes as fuel decreases but other dimensions do not change
   */
  void updateCOM(int dltaH)
  {
    com[Z] -= dltaH; 
  }
  int calcDltaH(int dltaM) {
    //TODO
    return 0;
  }
private:
  Shape s; //simple shape with known area, volume, other parameters, etc
  int massF; //just mass of fuel (mass of tank is included in dry_mass of rocket)
  int[] com; // x and y component of cof is constant if shape is symmetric/ z decreases as fuel in tank decreases, relative to P
  
}


class Rocket {
public:
  Rocket() {
    
  }
  /*
   * updatePARAM vs getPARAM (sometimes you just want to get the PARAM vs getting and updating it) this is done to avoid recalculating values that we already know
   *
   */
  double* updateCOM()
  {
    com[X] = dry_mass*length;
    com[Y] = dry_mass*width;
    com[Z] = dry_mass*height;

    for (int i = 0; i < tanks.length; i++)
    {
      for (int j = 0; j < DIM; j++)
      {
        com[j] += tanks[i].mass*tanks[i].cof[j];
      }
    }
    return com;
  }
  double* updateCOP()
  {
    cop = cov;
    return cop;
  }
  double* getI()
  {
    return I;
  }
  int updateMass(float dt)
  {
    tot_mass = dry_mass;
    for(Tank t: tanks)
    {
      tot_mass += t.updateMass(dt);
    }
    return tot_mass;
  }
  /*
   * return Coeffient of drag based on the geometry of the surface of each direction of the rocket
   */
  float* getCD() {
    return cd;
  }
  double* getCOM() {
    return com;
  }
  double* getCOP() {
    return cop;
  }
  float* getSA() {
    return sa;
  }
  int getMass(){
    return tot_mass;
  }
  int getP() {
    return P;
  }
  
private:
  //constants 

  const int height; // meters
  const int radius; // meters
  const double[] cov; // center of volume
  const float[] cd; 
  //Shape[] s, eventually we will allow for arbitrary simple shapes but use rectangle for now
  const int P; //gimbal base, all rotations are taken with respect to this point
  const int dry_mass; 
  const float[] sa;
  const double[] I;
  const float mass_vel_out; 
  // variables



  int tot_mass; // Newtons , tot_mass = dry_mass + m_tank1 + m_tank2 + ... + m_tankN
  double[] com; //center of mass
  double[] cop; //center of pressure
  
  Tank[] tanks; //arbitary number of fuel tanks (lander should only have two but future designs may have different numbers)

};


 /* 
void RK4(double[] initialVec, double[] vec, double (*diffEq)(float, double, double), float dt) {
    int[] n = int[DIM];
 
    float[] k1, k2, k3, k4, k5;
 
    for (int i = 0; i < DIM; i++)
    {
        // Apply Runge Kutta Formulas to find
        // next value of y
        k1[i] = //dt*diffEq(0, initialVec);
        k2[i] = dt*diffEq(.5*dt, initialVec + 0.5*k1[i], );
        k3[i] = dt*diffEq(.5*dt, initialVec + 0.5*k2[i]);
        k4[i] = dt*diffEq(dt, initialVec + k3[i]);
 
        vec[i] += (1.0/6.0)*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
    }
} */

class Control {
  public:
  Control () {

  }
  double[] getThrust() {
    return thrust;
  }
  private:
    double[] thrust;
}
/*
 * ORDER of calculating values
 * ALWAYS update F before M 
 * 
 *  updateFdrag
 *  updateFg
 *  updateFThrust
 * 
 *  updateNetF
 * 
 *  updateMdrag
 *  updateMg
 *  updateMThrust
 * 
 *  updateNetM
 *  updateAcc
 * 
 *  updateNetM
 *  updateAngAcc
 *  
 * 
 *  clearNet
 */
class Dynamics { //FINISH THE CONSTRUCTOR
public:
  Dynamics () {
    ang[DIM]; //global frame variables
    pos[DIM];
    vel[DIM];
    ang_vel[DIM];
    acc[DIM];
    ang_acc[DIM];
    inerM; //inertial force creating a moment


    //List of forces on the rocket
    double[] fDrag = double[DIM]; // just update forces so we don't have to make new memory have time the function is called
    double[] fThrust = double[DIM];
    double[] fG = double[DIM];

    //List of moment on the rocket;
    double[] mDrag = double[DIM]; // just update forces so we don't have to make new memory have time the function is called
    double[] mThrust = double[DIM];
    double[] mG = double[DIM];


    double[] netF = double[DIM]; //net force
    double[] netM = double[DIM]; //net moment
    Rocket* r = new Rocket();
    Environment* e = new Environment();
    Control* c = new Control();
  }
  void RK4(float dt)
  {
    int[] n = int[DIM];
  
    double[] k1Pos, k2Pos, k3Pos, k4Pos;
    double[] k1Ang, k2Ang, k3Ang, k4Ang;  

//initial force calculations
    Dynamics* state0 = new Dymnamics(*this); //deep clone current state

    state0->updateInerMom();
    state0->updateFdrag();
    state0->updateFg();
    state0->updateFThrust();
    state0->updateNetF({fDrag,fThrust,fG});
    state0->updateAcc();

    state0->updateMdrag();
    state0->updateMg();
    state0->updateMThrust();
    state0->updateNetM({mDrag,mThrust,mG});
    state0->updateAngAcc();

    for (int i = 0; i < DIM; i++) {
      k1Pos[i] = dt*dPosdt(0*dt, r.getVel()[i], r.getAcc()[i]); //dltaPos
      k1Ang[i] = dt*dPosdt(0*dt, r.getAngVel()[i], r.getAngAcc()[i]); //dltaAng
    }
    

    //this order is important we update position based off .5*dltaPos
    //then we update vel because its already applied from dPos/dt
    state0->updatePos(.5, k1Pos);
    state0->updateVel(.5*dt);
    state0->updateAng(.5, k1Ang);
    state0->updateAngVel(.5*dt);

    state0->r->updateMass(.5*dt);
    state0->r->updateCOM();
    state0->updateInerMom();
    state0->updateFdrag();
    state0->updateFg();
    state0->updateNetF({fDrag,fThrust,fG});
    state0->updateAcc();

    state0->updateMdrag();
    state0->updateMg();
    state0->updateMThrust();
    state0->updateNetM({mDrag,mThrust,mG});
    state0->updateAngAcc();
    state0->clearAcc();

    for (int i = 0; i < DIM; i++) {
      k2Pos[i] = dt*dPosdt(.5*dt, r.getVel()[i], r.getAcc()[i]);
      k2Ang[i] = dt*dPosdt(.5*dt, r.getAngVel()[i], r.getAngAcc()[i]);
    }

    state0->updatePos(.5, k2Pos);
    state0->updateVel(.5*dt);
    state0->updateAng(.5, k2Ang);
    state0->updateAngVel(.5*dt);

    state0->r->updateMass(.5*dt);
    state0->r->updateCOM();
    state0->updateInerMom();
    state0->updateFdrag();
    state0->updateFg();
    state0->updateNetF({fDrag,fThrust,fG});
    state0->updateAcc();
    
    state0->updateMdrag();
    state0->updateMg();
    state0->updateMThrust();
    state0->updateNetM({mDrag,mThrust,mG});
    state0->updateAngAcc();
    state0->clearAcc();

    for (int i = 0; i < DIM; i++) {
      k3Pos[i] = dt*dPosdt(.5*dt, r.getVel()[i], r.getAcc()[i]);
      k3Ang[i] = dt*dPosdt(.5*dt, r.getAngVel()[i], r.getAngAcc()[i]);
    }

    state0->updatePos(dt, k3Pos);
    state0->updateVel(dt);
    state0->updateAng(dt, k3Ang);
    state0->updateAngVel(dt);

    state0->r->updateMass(dt);
    state0->r->updateCOM();
    state0->updateInerMom();
    state0->updateFdrag();
    state0->updateFg();
    state0->updateNetF({fDrag,fThrust,fG});
    state0->updateAcc();

    state0->updateMdrag();
    state0->updateMg();
    state0->updateMThrust();
    state0->updateNetM({mDrag,mThrust,mG});
    state0->updateAngAcc();
    state0->clearAcc();

    for (int i = 0; i < DIM; i++) {
      k4Pos[i] = dt*dPosdt(dt, r.getVel()[i], r.getAcc()[i]);
      k4Ang[i] = dt*dPosdt(dt, r.getAngVel()[i], r.getAngAcc()[i]);
    }

    this->pos[i] += (1.0/6.0)*(k1Pos[i] + 2*k2Pos[i] + 2*k3Pos[i] + k4Pos[i]);
    this->ang[i] += (1.0/6.0)*(k1Ang[i] + 2*k2Ang[i] + 2*k3Ang[i] + k4Ang[i]);
  }
  void updateState() {
    this->r->updateCOM();
    this->updateInerMom();
    this->updateFdrag();
    this->updateFg();
    this->updateFThrust();
    this->updateNetF({fDrag,fThrust,fG});
    this->updateAcc();

    this->updateMdrag();
    this->updateMg();
    this->updateMThrust();
    this->updateNetM({mDrag,mThrust,mG});
    this->updateAngAcc();
  }
  void euler() {
    this->updateState();

    for (int i = 0; i < DIM; i++)
    {
      pos[i] += dPosdt(dt, vel[i], acc[i]);
      vel[i] += acc[i]*dt;
    }
    for (int i = 0; i < DIM; i++)
    {
      ang[i] += dPosdt(dt, ang_vel[i], ang_acc[i]);
      ang_vel[i] += ang_acc[i]*dt;
    }
    clearAcc();
  }
  double dPosdt(float dt, double vel, double acc)
  {
    return .5*acc*Math.pow(dt,2) + vel*dt;
  }
  /*
   * this may get replaced/moved in the future but this succesfully calculates the cross_p, assumes DIM = 3
   */ 
  void crossProduct(double dist[], double force[], double moment[])
  {
    moment[0] = dist[1] * force[2] - dist[2] * force[1];
    moment[1] = dist[2] * force[0] - dist[0] * force[2];
    moment[2] = dist[0] * force[1] - dist[1] * force[0];
  }

  /*
   * @param double[] air, model air resistance as a vector 
   * vel is the velocity of the rocket, need to find relative wind vector
   * cd is coffiecient of drag of the rocket in each vector direction / it is proportional to surface area of that part of the rocket
   * TODO account for air resistance that flow around the rocket
   */
  void updateFdrag()
  {
    double P = getP(pos[Z]); 
    for (int i = 0; i < DIM; i++)
    {
      fDrag[i] = .5*Math.pow((-vel[i] + e->air[i]),2)*P*r->getCD()[i]*r->getSA()[i];
    } 
  }
  /*
   * the moment created from the force of drag
   */
  void updateMdrag()
  {
    double dist = double[DIM];
    for (int i = 0; i < DIM; i++) {
      dist[i] = r->getCOP()[i]-(r->getP()[i]);
    }
    crossProduct(dist, fDrag, mDrag);
  }



  void updateFg()
  {
    fG[Z] = getG(pos[Z])*r->getMass(); 
  }
  void updateMg()
  {
    double dist = double[DIM];
    for (int i = 0; i < DIM; i++) {
      dist[i] = r->getCOM()[i]-r->getP()[i];
    }
    crossProduct(dist, fG, mG);
  }

   /*
   * given by the control block
   */
  void updateFThrust()
  {
    for (int i = 0; i < DIM; i++)
    {
      fThrust[i] = c->thrust[i];
    }
  }
  void updateMThrust()
  {
    double dist = double[DIM];
    for (int i = 0; i < DIM; i++) {
      dist[i] = r->getCOM()[i]-r->getP()[i];
    }
    crossProduct(dist, fThrust, mThrust);
  }


  void updateInerMom()
  {
    double dist = double[DIM];
    for (int i = 0; i < DIM; i++) {
      dist[i] = r->getCOM()[i]-r->getP()[i];
    }
    crossProduct(dist, netF, inerM);
  }

  void updateNetF(double[][] forces)
  {
    for (int i = 0; i < forces.length; i++) {
      for (int j = 0; j < DIM; j++) {
        fNet[j] += forces[i][j]; 
      }
    }

  }
  void updateNetM(double[][] moments)
  {
    for (int i = 0; i < moments.length; i++) {
      for (int j = 0; j < DIM; j++) {
        mNet[j] += moments[i][j]; 
      }
    }
  } 

  void clearNet()
  {
    for (int i = 0; i < DIM; i++){
      fNet[i] = 0;
      mNet[i] = 0;
    }
  }
  void updateAcc()
  {
    for (int i = 0; i < DIM; i++){
      acc[i] = netF[i]/r.getMass();
    }
  }
  void updateVel(float dt) {
    for (int i = 0; i < DIM; i++)
    {
      vel[i] += dt*acc[i];
    }
  }
  void updatePos(float dt) {
    for (int i = 0; i < DIM; i++)
    {
      pos[i] += dt*vel[i];
    }
  }
  void updateVel(float dt, double[] acc) {
    for (int i = 0; i < DIM; i++)
    {
      vel[i] += dt*acc[i];
    }
  }
  void updatePos(float dt, double[] vel) {
    for (int i = 0; i < DIM; i++)
    {
      pos[i] += dt*vel[i];
    }
  }

  void updateAngAcc()
  {
    
    for (int i = 0; i < DIM; i++){
      ang_acc[i] = (netM[i]- inerM)/r.getI()[i];
    }
  }
  void updateAngVel(float dt) {
    for (int i = 0; i < DIM; i++)
    {
      ang_vel[i] += dt*ang_acc[i];
    }
  }
  void updateAng(float dt) {
    for (int i = 0; i < DIM; i++)
    {
      ang[i] += dt*ang_vel[i];
    }
  }


  
 
  
private:
  double[] ang; //global frame variables
  double[] pos;
  double[] vel;
  double[] ang_vel;
  double[] acc;
  double[] ang_acc;
  double[] inerM; //inertial force creating a moment


  //List of forces on the rocket
  double[] fDrag = double[DIM]; // just update forces so we don't have to make new memory have time the function is called
  double[] fThrust = double[DIM];
  double[] fG = double[DIM];

  //List of moment on the rocket;
  double[] mDrag = double[DIM]; // just update forces so we don't have to make new memory have time the function is called
  double[] mThrust = double[DIM];
  double[] mG = double[DIM];


  double[] netF = double[DIM]; //net force
  double[] netM = double[DIM]; //net moment
  Rocket* r = new Rocket();
  Environment* e = new Environment();
  Control* c = new Control();
};

class Environment {
public:

  const double GM = 3.986004418*Math.pow(10,14); 
  const double R = 5927399.424; //radius of Earth in San Diego (meters)
  double[] air = double[DIM];

  double getAir()
  {
    return air;
  }
  double getG(double h)
  {
    return GM/Math.pow((R + h),2);
  }
  /*
   * density (kg/m) of atmosphere as a function of height (meters)
   * https://www.grc.nasa.gov/www/k-12/rocket/atmosmet.html
   */
  double getP(int h)
  {
    if (h < 11000) {
      return 101.29*Math.pow(((15.04-.00649*h + 273.1)/288.08),5.256);
    }
    else if (h > 25000) {
      return 2.488*Math.pow(((-131.21+.00299*h + 273.1)/216.6),-11.388);    }
    else {
      return 22.65*Math.pow(Math.e,(1.73-.000157*h));
    }
  }
  
  Environment() {

  }
};

int main() {
  float dt = .1f; //.1 seconds per physics tick
  float timeElapsed = 0;
  float timeFinal = 1800; //(s) 30 minutes

  Rocket* r = new Rocket();
  Environment* e = new Environment();
  double[] ang; //global frame variables
  double[] pos;
  double[] vel;
  double[] ang_vel;
  double[] acc;
  double[] ang_acc;
  double[] inerM; //inertial force creating a moment


  Dynamics sym = new Dynamics();


  while (timeElapsed < timeFinal) {
    sym.euler(dt);

    timeElapsed += dt;
  }
  
}