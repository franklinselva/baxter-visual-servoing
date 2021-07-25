#include <visp/vpFeaturePoint.h>
#include <ecn_baxter_vs/baxter_arm.h>
#include <visp/vpSubMatrix.h>
#include <visp/vpSubColVector.h>
//#include <ecn_common/vpQuadProg.h>
#include <ecn_common/visp_utils.h>
#include <math.h>

using namespace std;

int main(int argc, char** argv)
{
    //BaxterArm arm(argc, argv);    // defaults to simulation with right arm
    BaxterArm arm(argc, argv, false, "left");   // real robot with left arm


    vpColVector q = arm.init(), q_d = arm.init(), q_smax = arm.init(), q_smin = arm.init();

    vpColVector qmin = arm.jointMin(),
                qmax = arm.jointMax();
	
	q_d[3] = (qmin[3]+qmax[3])/2;
	q_d[4] = (qmin[4]+qmax[4])/2;
	q_d[5] = (qmin[5]+qmax[5])/2;
	q_d[6] = (qmin[6]+qmax[6])/2;
	
    // define a simple 2D point feature and its desired value
    vpFeaturePoint p,pd;
    pd.set_xyZ(0,0,1);

    // the error
    vpColVector e(3), e_tilde(10);
    double x, y, area;
    // desired area
    const double area_d = arm.area_d();

    // loop variables
    vpColVector qdot(7), Larea(6);
    vpMatrix L(3, 6), Js(3,7), H(10, 10), J_tilde(10, 7), I, hDiag(7, 7);
	auto A = 0., B = 0., C = 1.;
	
	H.eye(10);
	I.eye(7);
	
    while(arm.ok())
        {
            cout << "-------------" << endl;

            // get point features
            x = arm.x();
            y = arm.y();

			auto lambda = arm.lambda();
			auto rho = arm.rho();
			
			q_smax = qmax - rho*(qmax - qmin);
			q_smin = qmin + rho*(qmax - qmin);
			
            area = arm.area();
            p.set_xyZ(x,y, 1);
            std::cout << "x: " << x << ", y: " << y << ", area: " << area << '\n';

            // update error vector e
            ecn::putAt(e,p.error(pd),0);
            e[2] = area - area_d;

            // update interaction matrix L
            Larea[0] = -area*A;
            Larea[1] = -area*B;
            Larea[2] = area*(3/1 - C);
            Larea[3] = 3*area*y;
            Larea[4] = -3*area*x;
            Larea[5] = 0;
            
			ecn::putAt(L, p.interaction(vpBasicFeature::FEATURE_ALL), 0, 0);
			ecn::putAt(L, Larea.t(), 2, 0);
			
            // compute feature Jacobian from L and cameraJacobian
			Js = L * arm.cameraJacobian();
						
            // build H matrix (2nd section) using arm.rho()
            q = arm.jointPosition();

			//Basic Visual Servoing
			qdot = -lambda * Js.pseudoInverse() * e;

			// H Matrix			
			for (int i = 0; i <= 6; i++)
			{
				if (q[i] >= q_smax[i])
				{
					double h = ecn::weight(q[i], q_smax[i], qmax[i]);
					H[i+3][i+3] = h;
				}
				else
				{
					double h = ecn::weight(q[i], q_smin[i], qmin[i]);
					H[i+3][i+3] = h;
				}
			}

			// Jacobian
			
			ecn::putAt(J_tilde, Js, 0, 0);
			ecn::putAt(J_tilde, I, 3, 0);
			
			//error matrix
			ecn::putAt(e_tilde, e, 0);
			ecn::putAt(e_tilde, q - q_d, 3);
			
			//qdot = (H*J_tilde*qdot + lambda * H * e_tilde)^2;	
			qdot = -lambda * (H * J_tilde).pseudoInverse() * H * e_tilde;
			
            // send this command to the robot
            arm.setJointVelocity(qdot);

            // display current joint positions and VS error
            arm.plot(e);
    }
}
