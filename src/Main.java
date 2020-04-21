import java.awt.Frame;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.ArrayList;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.GLAutoDrawable;
import com.jogamp.opengl.GLCapabilities;
import com.jogamp.opengl.GLEventListener;
import com.jogamp.opengl.GLProfile;
import com.jogamp.opengl.awt.GLCanvas;
import com.jogamp.opengl.glu.GLU;
import com.jogamp.opengl.util.Animator;


public class Main extends GLCanvas implements GLEventListener, MouseWheelListener {
	
			static double timeSync = 0.01646666666;
			static double timeAcceleration = 0.0000035;
	
			// Minhas constantes.
			static double currentTime = 0;
			static double massaDaTerra = 5.98*1E24;
			static double raioDaTerra = 6.38*1E6;
			static double deltaTime = timeAcceleration*timeSync;
			static double scale = 100;
			static double totalTime = 1;
			static double gravitationalConstant = 6.67384 * 1E-11;
			// Fim das minhas constantes.
			
	 
	private static final int height = 1000;
	private static final int width = height;
	
	
	//static Particle myParticleDemo("Name", Mass, Radius, PositionX, PositionY, PositionZ, VelocityX, VelocityY, VelocityX, SpinX, SpinY, SpinZ) 
	
	
	static Particle myParticle1 = new Particle("0",1,1,10,0,-40,0,0,0,0,0,0);

	
	static ArrayList<Particle> myParticles = new ArrayList<Particle>();
	
	public static void main(String[] args) {
		
		
		
		for (double i = -2 ; i < 3 ; i++) 	{
			for (double j = -2 ; j <3 ; j++) {
				Particle myParticle1 = new Particle("0",1,1,i*15,0,j*15,0,0,0,0,0,0);
				myParticles.add(myParticle1);
			}
		}
		
		
		
		GLProfile glp = GLProfile.getDefault();
		GLCapabilities capabilities = new GLCapabilities(glp);
        GLCanvas canvas = new GLCanvas(capabilities);
        
        
        Frame frame = new Frame("AWT Window Test");
        frame.setSize(height, width);
        frame.add(canvas);
        frame.setVisible(true);

        
        frame.addWindowListener(new WindowAdapter() {
            public void windowClosing(WindowEvent e) {
                System.exit(0);
            }});
        
        canvas.addGLEventListener(new Main());
        canvas.addMouseWheelListener(new Main());
        Animator animator = new Animator(canvas);
        animator.start();
	}
	
	

		private void update() {
			//System.out.println(currentTime);
			double tempX = 0;
			double tempY = 0;
			double tempZ = 0;
			double p = 0;
			double temperature = 0;
			
			
			
			for (int i = 0; i < myParticles.size(); i++){
			executeAlgorithm(myParticles, deltaTime, totalTime);
			
			
			for (int j = 0; j < myParticles.size() ; j++) {
				p += Math.sqrt(Math.pow(myParticles.get(j).getvX(),2)+Math.pow(myParticles.get(j).getvZ(),2));
			}
			
			for (int j = 0; j < myParticles.size() ; j++) {
				temperature += ((1 / (myParticles.size() - 1)*2) * myParticles.get(j).getMass()*(Math.pow(myParticles.get(j).getvX(),2)+Math.pow(myParticles.get(j).getvZ(),2)));
			}
			
			System.out.println(temperature);
			
			}
			
			for (int i = 0; i < myParticles.size(); i++){
				if(myParticles.get(i).getxX() > scale) {
					myParticles.get(i).setxX(myParticles.get(i).getxX()-(scale*2));
				}
				
				if(myParticles.get(i).getxX() < -scale) {
					myParticles.get(i).setxX(myParticles.get(i).getxX()+(scale*2));
				}
				
				if(myParticles.get(i).getxZ() > scale) {
					myParticles.get(i).setxZ(myParticles.get(i).getxZ()-(scale*2));
				}
				
				if(myParticles.get(i).getxZ() < -scale) {
					myParticles.get(i).setxZ(myParticles.get(i).getxZ()+(scale*2));
				}
				
				
				}
			
			
			 /*
			for (int i = 0 ; i < myParticles.size() ; i++) {
				tempX = tempX + myParticles.get(i).getMass()*myParticles.get(i).getvX();
				tempY = tempY + myParticles.get(i).getMass()*myParticles.get(i).getvY();
				tempZ = tempZ + myParticles.get(i).getMass()*myParticles.get(i).getvZ();
			}
			
			/*for (int i = 0 ; i < myParticles.size() ; i++) {
				temp = temp + myParticles.get(i).getMass()*(myParticles.get(i).getxX()*myParticles.get(i).getvZ()-myParticles.get(i).getxZ()*myParticles.get(i).getvX());
			}
			
			
			
			
			System.out.println(tempX + "\t\t" + tempY + "\t\t" + tempZ);
			*/
			
			
		}
		
			
			
		
			private void drawCircle (GL2 gl, double bodyRadius, double scale, double positionX, double positionY, double screenHeight , double screenWidth) {
			float red = (float)Math.abs(Math.sin(bodyRadius));
			float green = (float)Math.abs(Math.cos(bodyRadius));
			float blue = (float)Math.abs(Math.tan(bodyRadius));
			gl.glBegin(GL.GL_POINTS);
			gl.glColor3f(red, green, blue);
			bodyRadius = bodyRadius/scale;
		    for(int i =0; i <= 300; i++){
		    double angle = 2 * Math.PI * i / 300;
		    double x = Math.cos(angle);
		    double y = Math.sin(angle);
		    for (double j = 1 ; j >= 0.1 ; j = j - 0.02) {
		    gl.glVertex2d(positionX+((x*0.001*width)*j)*bodyRadius,positionY+((y*0.001*height)*j)*bodyRadius);
		    }
		    }
		}



			private void render(GLAutoDrawable drawable) {
				
			    GL2 gl =  drawable.getGL().getGL2();
			    
			    for (int i = 0 ; i < myParticles.size() ; i ++ ) {
				    gl = drawable.getGL().getGL2();
				    gl.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT);
				    drawCircle (gl, myParticles.get(i).getRadius(),scale,myParticles.get(i).getxX()/scale, myParticles.get(i).getxZ()/scale, height, width);

			    }
			    gl.glEnd();

			}
	
	public void display(GLAutoDrawable drawable) {
		 	update();
		    render(drawable);
		    
	}
	
	

	
	
	public void init(GLAutoDrawable drawable) {
	    // put your OpenGL initialization code here
		GLU glu = new GLU();
		 drawable.getGL().setSwapInterval(1);
		 
	}

	public void dispose(GLAutoDrawable drawable) {
	    // put your cleanup code here
	}

	public void reshape(GLAutoDrawable drawable, int x, int y, int width, int height) {
		GL gl = drawable.getGL();
        gl.glViewport(0, 0, width, height);
	}
	
	private static void executeAlgorithm(ArrayList<Particle> myParticles, double deltaTime, double totalTime ) {
			ArrayList<ArrayList<Double>> netForces = new ArrayList<ArrayList<Double>>();
			//double distance = Math.sqrt((Math.pow(myParticle1.getxX()-myParticle2.getxX(),2)+Math.pow(myParticle1.getxY()-myParticle2.getxY(),2)+Math.pow(myParticle1.getxZ()-myParticle2.getxZ(),2)));
			//double velocityModule = Math.sqrt((Math.pow(myParticle1.getvX()-myParticle2.getvX(),2)+Math.pow(myParticle1.getvY()-myParticle2.getvY(),2)+Math.pow(myParticle1.getvZ()-myParticle2.getvZ(),2)));
			
			
			//returnList.add(currentTime);
			//returnList.add(myParticle.getxX());
			//returnList.add(myParticle.getxY());
			//returnList.add(myParticle.getxZ());
			//returnList.add(distance-myPlanet.getRadius()-myParticle.getRadius());
			//returnList.add(((Math.pow(velocityModule,2)*myParticle.getMass())/2)+(9.8*myParticle.getMass()*(distance-myPlanet.getRadius())));
			//returnList.add(Math.sqrt(Math.pow(myParticle.getxX()-0,2)+Math.pow(myParticle.getxY()-0,2)+Math.pow(myParticle.getxZ()-myPlanet.getRadius()-myParticle.getRadius(),2)));
			//returnList.add(velocityModule);
			
			
			for (int i = 0; i < myParticles.size() ; i++) {
			netForces.add(applyForces(i,myParticles,deltaTime));
			}
			
			
			//System.out.println(netForces.get(0).get(0));;
			
			integrateForcesEulerCromer(netForces,myParticles,myParticles);
			
			
			//distance = Math.sqrt((Math.pow(myParticle1.getxX()-myParticle2.getxX(),2)+Math.pow(myParticle1.getxY()-myParticle2.getxY(),2)+Math.pow(myParticle1.getxZ()-myParticle2.getxZ(),2)));
			currentTime += deltaTime;

		
	}

		private static ArrayList<Double> applyForces(int i, ArrayList<Particle> myParticles, double deltaTime){
			ArrayList<ArrayList<Double>> tempList = new ArrayList<ArrayList<Double>>();
			ArrayList<Double> returnList = new ArrayList<Double>();

			double xNetForce = 0;
			double yNetForce = 0;
			double zNetForce = 0;


			for (int j = 0; j < myParticles.size(); j++) {
				if (i!=j) {
			//tempList.add(forceGravitational(myParticles.get(i), myParticles.get(j)));
			tempList.add(forceLennardJones(myParticles.get(i), myParticles.get(j)));
			//tempList.add(forceDrag(myParticle,myPlanet,1,0.1));
			//tempList.add(forceMagnus(myParticle, myPlanet));
			//tempList.add(forceSpring(myParticle1, myParticle2, 0, 0 , 0, 1));
			//tempList.add(forcePendulum(myParticle1, myParticle2, 0, 0, 10));
				}
			}
			
			for (ArrayList<Double> j: tempList){
				
			xNetForce = xNetForce + j.get(0);
			yNetForce = yNetForce + j.get(1);
			zNetForce = zNetForce + j.get(2);
			
			}
			
			
			returnList.add(xNetForce);
			returnList.add(yNetForce);
			returnList.add(zNetForce);

			return returnList;
		}
		
		 private static ArrayList<Particle> integrateForcesEuler(ArrayList<ArrayList<Double>> netForces, ArrayList<Particle> myParticles, ArrayList<Particle> cloneList) {
			 	ArrayList<Particle> returnList = new ArrayList<Particle>();

			 	for (int i = 0 ; i < myParticles.size() ; i++) {
			 	
			 	myParticles.get(i).setvX(myParticles.get(i).getvX() + (deltaTime*(netForces.get(i).get(0)/myParticles.get(i).getMass())));
				myParticles.get(i).setvY(myParticles.get(i).getvY() + (deltaTime*(netForces.get(i).get(1)/myParticles.get(i).getMass())));
				myParticles.get(i).setvZ(myParticles.get(i).getvZ() + (deltaTime*(netForces.get(i).get(2)/myParticles.get(i).getMass())));

				myParticles.get(i).setxX(myParticles.get(i).getxX() +  deltaTime*myParticles.get(i).getvX());
				myParticles.get(i).setxY(myParticles.get(i).getxY() +  deltaTime*myParticles.get(i).getvY());
				myParticles.get(i).setxZ(myParticles.get(i).getxZ() +  deltaTime*myParticles.get(i).getvZ());
				
				
			 	}
			 	
				return returnList;
		}
		 
		 
		 private static ArrayList<Particle> integrateForcesEulerCromer(ArrayList<ArrayList<Double>> netForces, ArrayList<Particle> myParticles, ArrayList<Particle> cloneList) {
			 	ArrayList<Particle> returnList = new ArrayList<Particle>();
			 	
			 	double vX, vY, vZ = 0;

			 	for (int i = 0 ; i < myParticles.size() ; i++) {
			 		
			 		vX = myParticles.get(i).getvX() + (deltaTime*(netForces.get(i).get(0)/myParticles.get(i).getMass()));
			 		vY = myParticles.get(i).getvY() + (deltaTime*(netForces.get(i).get(1)/myParticles.get(i).getMass()));
			 		vZ = myParticles.get(i).getvZ() + (deltaTime*(netForces.get(i).get(2)/myParticles.get(i).getMass()));
			 	
			 	myParticles.get(i).setvX(vX);
				myParticles.get(i).setvY(vY);
				myParticles.get(i).setvZ(vZ);

				myParticles.get(i).setxX(myParticles.get(i).getxX() +  deltaTime*vX);
				myParticles.get(i).setxY(myParticles.get(i).getxY() +  deltaTime*vY);
				myParticles.get(i).setxZ(myParticles.get(i).getxZ() +  deltaTime*vZ);
				
				
			 	}
			 	
				return returnList;
		}
		 
		 
		 
		 private static void integrateForcesEulerRichardson(ArrayList<ArrayList<Double>> netForces, ArrayList<Particle> myParticles, ArrayList<Particle> cloneList){
				double aX, aY, aZ = 0;
				double vX, vY, vZ = 0;
				double xX, xY, xZ = 0;
				double xNetForce = 0;
				double yNetForce = 0;
				double zNetForce = 0;
				
				
				for (int i = 0 ; i < myParticles.size() ; i++) {
					ArrayList<ArrayList<Double>> tempList = new ArrayList<ArrayList<Double>>();
					aX = netForces.get(i).get(0) / myParticles.get(i).getMass();
					aY = netForces.get(i).get(1) / myParticles.get(i).getMass();
					aZ = netForces.get(i).get(2) / myParticles.get(i).getMass();
					
					vX = myParticles.get(i).getvX() + 0.5*aX*deltaTime;
					vY = myParticles.get(i).getvY() + 0.5*aY*deltaTime;
					vZ = myParticles.get(i).getvZ() + 0.5*aZ*deltaTime;
					
					xX = myParticles.get(i).getxX() + 0.5*vX*deltaTime;
					xY = myParticles.get(i).getxY() + 0.5*vY*deltaTime;
					xZ = myParticles.get(i).getxZ() + 0.5*vZ*deltaTime;
					
					
					Particle tempParticle = new Particle(myParticles.get(i).getMass(),myParticles.get(i).getRadius(),xX,xY,xZ,vX,vY,vZ,myParticles.get(i).getSpinX(),myParticles.get(i).getSpinY(),myParticles.get(i).getSpinZ());
					
					for (int j = 0 ; j < myParticles.size(); j++) {
						if (i!=j){
							tempList.add(forceGravitational(tempParticle, myParticles.get(j)));
						}
					}
					
					System.out.println(tempList.get(0).get(0));
					
					for (ArrayList<Double> j: tempList){
						
						xNetForce = xNetForce + j.get(0);
						yNetForce = yNetForce + j.get(1);
						zNetForce = zNetForce + j.get(2);
						
						}
					
					aX = xNetForce / myParticles.get(i).getMass();
					aY = yNetForce / myParticles.get(i).getMass();
					aZ = zNetForce / myParticles.get(i).getMass();
					
					myParticles.get(i).setvX(myParticles.get(i).getvX()+aX*deltaTime);
					myParticles.get(i).setvY(myParticles.get(i).getvY()+aY*deltaTime);
					myParticles.get(i).setvZ(myParticles.get(i).getvZ()+aZ*deltaTime);
					
					myParticles.get(i).setxX(myParticles.get(i).getxX() + vX*deltaTime);
					myParticles.get(i).setxY(myParticles.get(i).getxY() + vY*deltaTime);
					myParticles.get(i).setxZ(myParticles.get(i).getxZ() + vZ*deltaTime);
					
					
				}
				
				
				
				
				
				
				
				
			}



		private static ArrayList<Double> forceSpring(Particle myParticle1, Particle myParticle2, double springX, double springY, double springZ, double k) {
			 ArrayList<Double> returnList = new ArrayList<Double>();
			 double forceX, forceY, forceZ;
			 double distance = Math.sqrt((Math.pow(myParticle1.getxX()-springX,2)+Math.pow(myParticle1.getxY()-springY,2)+Math.pow(myParticle1.getxZ()-springZ,2)));
			 double x = Math.sqrt(Math.pow(myParticle1.getxX()-springX,2)+Math.pow(myParticle1.getxY()-springY,2)+Math.pow(myParticle1.getxZ()-springZ,2));
			 double temp = -k*x;
			 
			 
			 
			 if (distance != 0) {
				 	forceX = (myParticle1.getxX()-springX)/distance;
					forceY = (myParticle1.getxY()-springY)/distance;
					forceZ = (myParticle1.getxZ()-springZ)/distance;
				 	} else {
				 		forceX = 0;
				 		forceY = 0;
				 		forceZ = 0;
				 	}
				 	
			 	forceX = forceX*temp;
				forceY = forceY*temp;
				forceZ = forceZ*temp;
			 
				returnList.add(forceX);
				returnList.add(forceY);
				returnList.add(forceZ);
			 
			 return returnList;
			 
			 
		 }
		 
		 private static ArrayList<Double> forcePendulum(Particle myParticle1, Particle myParticle2, double pendulumX, double pendulumY, double pendulumZ){
			 ArrayList<Double> returnList = new ArrayList<Double>();
			 double forceX, forceY, forceZ;
			 double distance = Math.sqrt((Math.pow(myParticle1.getxX()-myParticle2.getxX(),2)+Math.pow(myParticle1.getxY()-myParticle2.getxY(),2)+Math.pow(myParticle1.getxZ()-myParticle2.getxZ(),2)));
			 double gravity =(-gravitationalConstant * myParticle2.getMass()*myParticle1.getMass() )/ Math.pow(distance,2);
			 double velocityModule = Math.sqrt((Math.pow(myParticle1.getvX()-myParticle2.getvX(),2)+Math.pow(myParticle1.getvY()-myParticle2.getvY(),2)+Math.pow(myParticle1.getvZ()-myParticle2.getvZ(),2)));
			 double temp;
			 double stringSize = 10;
			 double cosineTheta = ((-1)*(myParticle1.getxZ()-pendulumZ))/stringSize;
			 
			 
			 
			 temp = -((myParticle1.getMass()*Math.pow(velocityModule,2))/stringSize) + (cosineTheta * myParticle1.getMass() * gravity);
			 
			 if (distance != 0) {
				 	forceX = (myParticle1.getxX()-pendulumX)/stringSize;
					forceY = (myParticle1.getxY()-pendulumY)/stringSize;
					forceZ = (myParticle1.getxZ()-pendulumZ)/stringSize;
				 	} else {
				 		forceX = 0;
				 		forceY = 0;
				 		forceZ = 0;
				 	}
				 	
			 	forceX = forceX*temp;
				forceY = forceY*temp;
				forceZ = forceZ*temp;
			 
				returnList.add(forceX);
				returnList.add(forceY);
				returnList.add(forceZ);
			 
				return returnList;
		 }
		
		 private static ArrayList<Double> forceGravitational(Particle myParticle1, Particle myParticle2){
				ArrayList<Double> returnList = new ArrayList<Double>();
				double temp;
				double forceX, forceY, forceZ, distance;
				distance = Math.sqrt((Math.pow(myParticle1.getxX()-myParticle2.getxX(),2)+Math.pow(myParticle1.getxY()-myParticle2.getxY(),2)+Math.pow(myParticle1.getxZ()-myParticle2.getxZ(),2)));

				temp =(-gravitationalConstant * myParticle2.getMass()*myParticle1.getMass() )/ Math.pow(distance,2);
				//System.out.println((distance-myPlanet.getRadius()));


				forceX = (myParticle1.getxX()-myParticle2.getxX())/distance;
				forceY = (myParticle1.getxY()-myParticle2.getxY())/distance;
				forceZ = (myParticle1.getxZ()-myParticle2.getxZ())/distance;
				forceX = forceX*temp;
				forceY = forceY*temp;
				forceZ = forceZ*temp;


				returnList.add(forceX);
				returnList.add(forceY);
				returnList.add(forceZ);
				
				return returnList;
			}
		 
		 private static ArrayList<Double> forceLennardJones(Particle myParticle1, Particle myParticle2){
				ArrayList<Double> returnList = new ArrayList<Double>();
				double temp;
				double forceX, forceY, forceZ, distance, distanceX, distanceY, distanceZ,x,y,z;
				distanceX = (myParticle1.getxX()-myParticle2.getxX());
				distanceY = (myParticle1.getxY()-myParticle2.getxY());
				distanceZ = (myParticle1.getxZ()-myParticle2.getxZ());
				x = 0;
				y=0;
				z=0;
				
				if (distanceX > scale) {
					x = 1;
				}
				
				if (distanceY > scale) {
					z = 1;
				}
				
				if (distanceZ > scale) {
					z = 1;
				}

				distanceX = pbcSeparation(distanceX,scale);
				distanceY = pbcSeparation(distanceY,scale);
				distanceZ = pbcSeparation(distanceZ,scale);
				
				
				distance = Math.sqrt((Math.pow(distanceX,2)+Math.pow(distanceY,2)+Math.pow(distanceZ,2)));
				temp =  (24/distance)* (2* Math.pow((scale/distance),12)- Math.pow((scale/distance),6));


				forceX = (distanceX)/distance;
				forceY = (distanceY)/distance;
				forceZ = (distanceZ)/distance;
				forceX = forceX*temp;
				forceY = forceY*temp;
				forceZ = forceZ*temp;


				
				returnList.add(forceX);
				
				
				returnList.add(forceY);
				
				
				returnList.add(forceZ);
				
				return returnList;
			}


		 	static private double pbcSeparation (double ds , double L) {
			 if ( ds > L) {
			 ds -= L*2 ;
			 } else if ( ds < -L) {
			 ds += L*2 ;
			 }
			 return ds ;
			 }
		 
		 
		@Override
		  public void mouseWheelMoved(MouseWheelEvent e) {
		       int notches = e.getWheelRotation();
		       double wheelFactor = 0.90;
		       if (notches < 0) {
		           scale = scale*0.5;
		       } else {
		           scale = scale*2;
		       }
		      
		  
		}
		 
}