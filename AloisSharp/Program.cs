using System;
using System.Collections.Generic;

namespace AloisSharp
{
    class Program
    {

        //Class load case
        public class LoadCase
        {

            public string Type; //Load type

            public List<unifLoad> ul;

            public List<nodeLoad> nl;

            public List<double>[] N;
            public List<double>[] M;
            public List<double>[] V;
            public List<double>[] Du;
            public List<double>[] Dv;


            public LoadCase(string Type)
            {
                this.Type = Type;

                ul = new List<unifLoad>();
                nl = new List<nodeLoad>();

            }

            public void Init(int nel)
            {
                N = new List<double>[nel];
                M = new List<double>[nel];
                V = new List<double>[nel];
                Du = new List<double>[nel];
                Dv = new List<double>[nel];


                for (int i = 0; i < nel; i++)
                {
                    N[i] = new List<double>();
                    M[i] = new List<double>();
                    V[i] = new List<double>();
                    Du[i] = new List<double>();
                    Dv[i] = new List<double>();

                }

            }


            public void addUnifLoad(unifLoad ul)
            {
                this.ul.Add(ul);
            }

            public void addNodeLoad(nodeLoad nl)
            {
                this.nl.Add(nl);
            }


            public double[][] genFle(List<El2d> Elements)
            {

                double[][] fle = new double[Elements.Count][];

                for (int i = 0; i < Elements.Count; i++)
                    fle[i] = new double[6];

                for (int i = 0; i < ul.Count; i++)
                {

                    fle[ul[i].nel][0] += Elements[ul[i].nel].getFle(ul[i])[0];
                    fle[ul[i].nel][1] += Elements[ul[i].nel].getFle(ul[i])[1];
                    fle[ul[i].nel][2] += Elements[ul[i].nel].getFle(ul[i])[2];
                    fle[ul[i].nel][3] += Elements[ul[i].nel].getFle(ul[i])[3];
                    fle[ul[i].nel][4] += Elements[ul[i].nel].getFle(ul[i])[4];
                    fle[ul[i].nel][5] += Elements[ul[i].nel].getFle(ul[i])[5];

                }

                return fle;

            }

        }

        //Class load combination
        public class LoadComb
        {
            public string loadComb; //Load type

            public List<List<double>> N;
            public List<List<double>> M;
            public List<List<double>> V;

            public LoadComb(string loadComb)
            {
                this.loadComb = loadComb;
            }

            public LoadComb()
            {

            }

        }

        //Class boundary condition
        public class BC
        {
            public int nn; //node number
            public int vx; //horizontal     1:fixed
            public int vy; //vertical       1:fixed
            public int rr; //rotation       1:fixed

            public BC(int nn, int vx, int vy, int rr)
            {
                this.nn = nn;
                this.vx = vx;
                this.vy = vy;
                this.rr = rr;
            }

            public BC()
            {
                vx = 0;
                vy = 0;
                rr = 0;
            }

        }

        //Class material
        public class Material
        {

            public double E;
            public double Gm;

            public Material(double E, double Gm)
            {
                this.E = E;
                this.Gm = Gm;
            }
        }

        //Class section
        public class Section
        {

            public double A;
            public double I;
            public double ks;

            public Section(double A, double I, double ks)
            {
                this.A = A;
                this.I = I;
                this.ks = ks;
            }

        }

        //Distribuited load
        public class unifLoad
        {
            public int nel;
            public double qx;
            public double qy;

            public unifLoad(int nel, double qx, double qy)
            {
                this.nel = nel;
                this.qx = qx;
                this.qy = qy;
            }

        }

        //Concentrated nodal load 
        public class nodeLoad
        {

            public int nn;
            public double fx, fy, m;

            public nodeLoad(int nn, double fx, double fy, double m)
            {
                this.nn = nn;
                this.fx = fx;
                this.fy = fy;
                this.m = m;
            }

            public nodeLoad()
            {
                fx = 0;
                fy = 0;
                m = 0;
            }


        }

        //Class beam element - Timoshenko 
        public class El2d
        {

            public int num;
            public double L;
            public double[] b;
            public double[] n;
            public double[][] K, G;
            public Material mat;
            public Section sec;
            public Node2d[] Nodes;

            public double[] fle;

            public string ElType;


            //Constructor
            public El2d(int num, string ElType, Node2d n1, Node2d n2, Material mat, Section sec)
            {

                this.num = num;

                this.ElType = ElType;

                Nodes = new Node2d[2];
                Nodes[0] = n1;
                Nodes[1] = n2;

                b = new double[2];
                b[0] = n2.x - n1.x;
                b[1] = n2.y - n1.y;

                L = Math.Sqrt(Math.Pow(b[0], 2) + Math.Pow(b[1], 2));

                n = new double[2];
                n[0] = b[0] / L;
                n[1] = b[1] / L;

                this.mat = mat;
                this.sec = sec;

                double E, Gm, A, I, ks, m;

                E = mat.E;
                Gm = mat.Gm;
                I = sec.I;
                A = sec.A;
                ks = sec.ks;

                m = 12 / Math.Pow(L, 2) * (E * I / (Gm * A * ks));

                K = new double[6][];

                //fle = new double[6];

                for (int i = 0; i < 6; i++)
                    K[i] = new double[6];

                if (ElType == "Beam2t")
                {  //Timoshenko beam element

                    K[0][0] = E / (1 + m) * A * (1 + m) / L;
                    K[1][1] = E / (1 + m) * 12 * I / Math.Pow(L, 3);
                    K[2][2] = E / (1 + m) * 4 * I * (1 + m / 4) / L;
                    K[3][3] = E / (1 + m) * A * (1 + m) / L;
                    K[4][4] = E / (1 + m) * 12 * I / Math.Pow(L, 3);
                    K[5][5] = E / (1 + m) * 4 * I * (1 + m / 4) / L;
                    K[0][1] = K[1][0] = 0;
                    K[0][2] = K[2][0] = 0;
                    K[0][3] = K[3][0] = E / (1 + m) * (-A * (1 + m) / L);
                    K[0][4] = K[4][0] = 0;
                    K[0][5] = K[5][0] = 0;
                    K[1][2] = K[2][1] = E / (1 + m) * 6 * I / Math.Pow(L, 2);
                    K[1][3] = K[3][1] = 0;
                    K[1][4] = K[4][1] = E / (1 + m) * (-12 * I / Math.Pow(L, 3));
                    K[1][5] = K[5][1] = E / (1 + m) * 6 * I / Math.Pow(L, 2);
                    K[2][3] = K[3][2] = 0;
                    K[2][4] = K[4][2] = E / (1 + m) * (-6 * I / Math.Pow(L, 2));
                    K[2][5] = K[5][2] = E / (1 + m) * 2 * I * (1 - m / 2) / L;
                    K[3][4] = K[4][3] = 0;
                    K[3][5] = K[5][4] = 0;
                    K[4][5] = K[5][4] = E / (1 + m) * (-6 * I / Math.Pow(L, 2));

                }
                else if (ElType == "Beam2t_hr")
                {   //Timoshenko beam element - right hinge

                    K[0][0] = E / (1 + m) * A * (1 + m) / L;
                    K[1][1] = 12 * E * I / (Math.Pow(L, 3) * (m + 1)) - 36 * E * I / (Math.Pow(L, 3) * (m + 4) * (m + 1));
                    K[2][2] = E * I * (m + 4) / (L * (m + 1)) - E * I * Math.Pow((m - 2), 2) / (L * (m + 4) * (m + 1));
                    K[3][3] = E / (1 + m) * A * (1 + m) / L;
                    K[4][4] = 12 * E * I / (Math.Pow(L, 3) * (m + 1)) - 36 * E * I / (Math.Pow(L, 3) * (m + 4) * (m + 1));
                    K[5][5] = 0;
                    K[0][1] = K[1][0] = 0;
                    K[0][2] = K[2][0] = 0;
                    K[0][3] = K[3][0] = E / (1 + m) * (-A * (1 + m) / L);
                    K[0][4] = K[4][0] = 0;
                    K[0][5] = K[5][0] = 0;
                    K[1][2] = K[2][1] = 6 * E * I / (Math.Pow(L, 2) * (m + 1)) + 6 * E * I * (m - 2) / (Math.Pow(L, 2) * (m + 4) * (m + 1));
                    K[1][3] = K[3][1] = 0;
                    K[1][4] = K[4][1] = -12 * E * I / (Math.Pow(L, 3) * (m + 1)) + 36 * E * I / (Math.Pow(L, 3) * (m + 4) * (m + 1));
                    K[1][5] = K[5][1] = 0;
                    K[2][3] = K[3][2] = 0;
                    K[2][4] = K[4][2] = -6 * E * I / (Math.Pow(L, 2) * (m + 1)) - 6 * E * I * (m - 2) / (Math.Pow(L, 2) * (m + 4) * (m + 1));
                    K[2][5] = K[5][2] = 0;
                    K[3][4] = K[4][3] = 0;
                    K[3][5] = K[5][3] = 0;
                    K[4][5] = K[5][4] = 0;

                }
                else if (ElType == "Beam2t_hl")
                {   //Timoshenko beam element - left hinge

                    K[0][0] = A * E / L;
                    K[1][1] = 12 * E * I / (Math.Pow(L, 3) * (m + 1)) - 36 * E * I / (Math.Pow(L, 3) * (m + 4) * (m + 1));
                    K[2][2] = 0;
                    K[3][3] = A * E / L;
                    K[4][4] = 12 * E * I / (Math.Pow(L, 3) * (m + 1)) - 36 * E * I / (Math.Pow(L, 3) * (m + 4) * (m + 1));
                    K[5][5] = E * I * (m + 4) / (L * (m + 1)) - E * I * Math.Pow((m - 2), 2) / (L * (m + 4) * (m + 1));
                    K[0][1] = K[1][0] = 0;
                    K[0][2] = K[2][0] = 0;
                    K[0][3] = K[3][0] = E / (1 + m) * (-A * (1 + m) / L);
                    K[0][4] = K[4][0] = 0;
                    K[0][5] = K[5][0] = 0;
                    K[1][2] = K[2][1] = 0;
                    K[1][3] = K[3][1] = 0;
                    K[1][4] = K[4][1] = -12 * E * I / (Math.Pow(L, 3) * (m + 1)) + 36 * E * I / (Math.Pow(L, 3) * (m + 4) * (m + 1));
                    K[1][5] = K[5][1] = 6 * E * I / (Math.Pow(L, 2) * (m + 1)) + 6 * E * I * (m - 2) / (Math.Pow(L, 2) * (m + 4) * (m + 1));
                    K[2][3] = K[3][2] = 0;
                    K[2][4] = K[4][2] = 0;
                    K[2][5] = K[5][2] = 0;
                    K[3][4] = K[4][3] = 0;
                    K[3][5] = K[5][3] = 0;
                    K[4][5] = K[5][4] = -6 * E * I / (Math.Pow(L, 2) * (m + 1)) - 6 * E * I * (m - 2) / (Math.Pow(L, 2) * (m + 4) * (m + 1));

                }
                else if (ElType == "Bar2d")
                {   //Truss element

                    K[0][0] = A * E / L;
                    K[1][1] = 0;
                    K[2][2] = 1; //Check Value!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    K[3][3] = A * E / L;
                    K[4][4] = 0;
                    K[5][5] = 1; //Check Value!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    K[0][1] = K[1][0] = 0;
                    K[0][2] = K[2][0] = 0;
                    K[0][3] = K[3][0] = E * (-A / L);
                    K[0][4] = K[4][0] = 0;
                    K[0][5] = K[5][0] = 0;
                    K[1][2] = K[2][1] = 0;
                    K[1][3] = K[3][1] = 0;
                    K[1][4] = K[4][1] = 0;
                    K[1][5] = K[5][1] = 0;
                    K[2][3] = K[3][2] = 0;
                    K[2][4] = K[4][2] = 0;
                    K[2][5] = K[5][2] = 0;
                    K[3][4] = K[4][3] = 0;
                    K[3][5] = K[5][3] = 0;
                    K[4][5] = K[5][4] = 0;

                }

                G = new double[6][];

                for (int i = 0; i < 6; i++)
                    G[i] = new double[6];

                G[0][0] = n[0];
                G[1][1] = n[0];
                G[2][2] = 1;
                G[3][3] = n[0];
                G[4][4] = n[0];
                G[5][5] = 1;
                G[0][1] = n[1];
                G[0][2] = 0;
                G[0][3] = 0;
                G[0][4] = 0;
                G[0][5] = 0;
                G[1][2] = 0;
                G[1][3] = 0;
                G[1][4] = 0;
                G[1][5] = 0;
                G[2][3] = 0;
                G[2][4] = 0;
                G[2][5] = 0;
                G[3][4] = n[1];
                G[3][5] = 0;
                G[4][5] = 0;
                G[1][0] = -n[1];
                G[2][0] = 0;
                G[3][0] = 0;
                G[4][0] = 0;
                G[5][0] = 0;
                G[2][1] = 0;
                G[3][1] = 0;
                G[4][1] = 0;
                G[5][1] = 0;
                G[3][2] = 0;
                G[4][2] = 0;
                G[5][2] = 0;
                G[4][3] = -n[1];
                G[5][3] = 0;
                G[5][4] = 0;

                K = Matrix.Multiply(Matrix.Multiply(Matrix.Traspose(G), K), G);

            }


            public double[] getFle(unifLoad ul)
            {

                fle = new double[6];

                if (ElType == "Beam2t")
                {  //Timoshenko beam element

                    fle[0] += L * ul.qx / 2;
                    fle[1] += L * ul.qy / 2;
                    fle[2] += L * ul.qy * L / 12;
                    fle[3] += L * ul.qx / 2;
                    fle[4] += L * ul.qy / 2;
                    fle[5] += -L * ul.qy * L / 12;

                }
                else if (ElType == "Beam2t_hr")
                {   //Timoshenko beam element - right hinge

                    fle[0] += L * ul.qx / 2;
                    fle[1] += 5.0 / 8.0 * L * ul.qy;
                    fle[2] += L * ul.qy * L / 8;
                    fle[3] += L * ul.qx / 2;
                    fle[4] += 3.0 / 8.0 * L * ul.qy;
                    fle[5] += 0;

                }
                else if (ElType == "Beam2t_hl")
                {   //Timoshenko beam element - left hinge

                    fle[0] += L * ul.qx / 2;
                    fle[1] += 3.0 / 8.0 * L * ul.qy;
                    fle[2] += 0;
                    fle[3] += L * ul.qx / 2;
                    fle[4] += 5.0 / 8.0 * L * ul.qy;
                    fle[5] += -L * ul.qy * L / 8;

                }
                else if (ElType == "Bar2d")
                {   //Truss element

                    fle[0] += L * ul.qx / 2;
                    fle[1] += 0;
                    fle[2] += 0;
                    fle[3] += L * ul.qx / 2;
                    fle[4] += 0;
                    fle[5] += 0;

                }

                fle = Matrix.Multiply(Matrix.Traspose(G), fle);

                return fle;
            }
        }

        //Class node in 2d
        public class Node2d
        {

            public int num;
            public double x;
            public double y;

            public int ndof1;
            public int ndof2;
            public int ndof3;

            public BC bc;

            //public List<nodeLoad> nl; 


            //Constructor
            public Node2d(double x, double y)
            {

                this.x = x;
                this.y = y;

                bc = new BC();
                //nl = new List<nodeLoad>();

            }

            //Constructor
            public Node2d(int num, double x, double y)
            {

                this.x = x;
                this.y = y;
                this.num = num;

                ndof1 = 3 * num - 2;
                ndof2 = 3 * num - 1;
                ndof3 = 3 * num - 0;

                bc = new BC();
                //nl = new List<nodeLoad>();

            }

            /*
            public void addNodeLoad(nodeLoad nl)
            {
                this.nl.Add(nl);
            }
            */

        }

        //Class structure 2d (frame)
        public class Str2d
        {

            public double[][] K;

            public List<Node2d> Nodes;
            public List<El2d> Elements;

            public List<List<int>> edof;

            public double[][] l;
            public double[][] u;
            public double[] disp;
            public double[][] Kbc;

            public double[] f;
            public double[] f0;

            int ndof;

            double[][] ed;

            double[] React; //React



            //Constructor
            public Str2d()
            {

                this.Nodes = new List<Node2d>();
                this.Elements = new List<El2d>();
                edof = new List<List<int>>();


            }


            //Extract element displacements from the global displacement 
            //vector according to the topology matrix edof
            public void extract()
            {

                int nie, n;
                nie = edof.Count;
                n = edof[0].Count;

                ed = new double[nie][];

                for (int i = 0; i < nie; i++)
                    ed[i] = new double[n];

                for (int i = 0; i < nie; i++)
                    for (int j = 0; j < n; j++)
                        ed[i][j] = disp[edof[i][j] - 1];

            }

            //Internal forces Timoshenko
            public void IntForces(LoadCase lc, int el, int nn)
            {

                string ElType = Elements[el].ElType;

                if (nn < 2) nn = 2;

                double EA, EI, GAK, alfa, L;
                double qx = 0.0, qy = 0.0;
                double[] n = new double[2];
                double[] b = new double[2];
                double[][] G = new double[6][];
                double[][] invC = new double[6][];

                EA = Elements[el].mat.E * Elements[el].sec.A;
                EI = Elements[el].mat.E * Elements[el].sec.I;
                GAK = Elements[el].mat.Gm * Elements[el].sec.A * Elements[el].sec.ks;
                alfa = EI / GAK;

                b[0] = Elements[el].Nodes[1].x - Elements[el].Nodes[0].x;
                b[1] = Elements[el].Nodes[1].y - Elements[el].Nodes[0].y;

                L = Math.Sqrt(Math.Pow(b[0], 2) + Math.Pow(b[1], 2));

                for (int i = 0; i < lc.ul.Count; i++)
                {
                    if (lc.ul[i].nel == el)
                    {
                        qx += lc.ul[i].qx;
                        qy += lc.ul[i].qy;
                    }
                }

                n[0] = b[0] / L;
                n[1] = b[1] / L;

                double[][] C = new double[6][];

                C[0] = new double[6] { 0, 0, 0, 1, 0, 0 };
                C[1] = new double[6] { 0, 0, 0, 0, 0, 1 };
                C[2] = new double[6] { 0, 6 * alfa, 0, 0, 1, 0 };
                C[3] = new double[6] { L, 0, 0, 1, 0, 0 };
                C[4] = new double[6] { 0, Math.Pow(L, 3), Math.Pow(L, 2), 0, L, 1 };
                C[5] = new double[6] { 0, 3 * (Math.Pow(L, 2) + 2 * alfa), 2 * L, 0, 1, 0 };

                invC = Matrix.inverse(C);

                G[0] = new double[6] { n[0], n[1], 0, 0, 0, 0 };
                G[1] = new double[6] { -n[1], n[0], 0, 0, 0, 0 };
                G[2] = new double[6] { 0, 0, 1, 0, 0, 0 };
                G[3] = new double[6] { 0, 0, 0, n[0], n[1], 0 };
                G[4] = new double[6] { 0, 0, 0, -n[1], n[0], 0 };
                G[5] = new double[6] { 0, 0, 0, 0, 0, 1 };


                double[] arr = new double[6] { 0, 0, 0, -qx * Math.Pow(L, 2) / (2 * EA), qy * Math.Pow(L, 4) / (24 * EI) - qy * Math.Pow(L, 2) / (2 * GAK), qy * Math.Pow(L, 3) / (6 * EI) };


                if (ElType == "Beam2t_hr")  //Timoshenko beam element - right hinge
                    arr = new double[6] { -qx * Math.Pow(L, 2) / (2 * EA), qy * Math.Pow(L, 4) / (24 * EI) - qy * Math.Pow(L, 2) / (2 * GAK), -qy * Math.Pow(L, 3) / (6 * EI), 0, 0, 0 };
                else if (ElType == "Beam2t_hl") //Timoshenko beam element - left hinge
                    arr = new double[6] { 0, 0, 0, -qx * Math.Pow(L, 2) / (2 * EA), qy * Math.Pow(L, 4) / (24 * EI) - qy * Math.Pow(L, 2) / (2 * GAK), qy * Math.Pow(L, 3) / (6 * EI) };



                double[] m = new double[6];

                m = Matrix.Multiply(G, ed[el]);

                m[0] = m[0] - arr[0];
                m[1] = m[1] - arr[1];
                m[2] = m[2] - arr[2];
                m[3] = m[3] - arr[3];
                m[4] = m[4] - arr[4];
                m[5] = m[5] - arr[5];

                double[] mm;
                mm = Matrix.Multiply(invC, m);

                double[] C2 = new double[2];
                C2[0] = mm[0];
                C2[1] = mm[3];

                double[] C4 = new double[4];
                C4[0] = mm[1];
                C4[1] = mm[2];
                C4[2] = mm[4];
                C4[3] = mm[5];

                double[] x = new double[nn];

                double ll = L / (nn - 1);

                for (int i = 0; i < nn; i++)
                    x[i] = i * ll;

                double[] u = new double[nn];
                double[] du = new double[nn];
                double[] teta = new double[nn];
                double[] dteta = new double[nn];
                double[] ddteta = new double[nn];

                for (int i = 0; i < nn; i++)
                    u[i] = x[i] * C2[0] + C2[1] - qx / (2 * EA) * Math.Pow(x[i], 2);

                for (int i = 0; i < nn; i++)
                    du[i] = C2[0] - qx * x[i] / (2 * EA);

                double[] v = new double[nn];
                double[] dv = new double[nn];

                //Timoshenko element without internal hinge
                if (ElType == "Beam2t")
                {

                    double v1, t1, v2, t2;
                    double a0, a1, a2, a3;
                    double ph, phs, phss;

                    ph = 12 * EI / (GAK * Math.Pow(L, 2));
                    phs = 1 / (1 + ph);
                    phss = EI / GAK;

                    v1 = m[1];
                    t1 = m[2];
                    v2 = m[4];
                    t2 = m[5];

                    a0 = v1;
                    a1 = phs * (-1 / L * ph * v1 + (1 + ph / 2) * t1 + 1 / L * ph * v2 - ph / 2 * t2);
                    a2 = phs * (-3 / Math.Pow(L, 2) * v1 - 1 / L * (2 + ph / 2) * t1 + 3 / Math.Pow(L, 2) * v2 - 1 / L * (1 - ph / 2) * t2);
                    a3 = phs * (2 * v1 / Math.Pow(L, 3) + t1 / Math.Pow(L, 2) - 2 * v2 / Math.Pow(L, 3) + t2 / Math.Pow(L, 2));

                    for (int i = 0; i < nn; i++)
                        v[i] = Math.Pow(x[i], 3) * a3 + Math.Pow(x[i], 2) * a2 + x[i] * a1 + a0 + qy / (24 * EI) * Math.Pow(x[i], 4) - qy / (2 * GAK) * Math.Pow(x[i], 2);

                    for (int i = 0; i < nn; i++)
                        dv[i] = 3 * Math.Pow(x[i], 2) * a3 + 2 * x[i] * a2 + a1 + qy * Math.Pow(x[i], 3) / (6 * EI) - qy * x[i] / GAK;

                    for (int i = 0; i < nn; i++)
                        teta[i] = 3 * (Math.Pow(x[i], 2) + 2 * alfa) * a3 + 2 * x[i] * a2 + a1 + qy * Math.Pow(x[i], 3) / (6 * EI);

                    for (int i = 0; i < nn; i++)
                        dteta[i] = 6 * x[i] * a3 + 2 * a2 + qy * Math.Pow(x[i], 2) / (2 * EI);

                    for (int i = 0; i < nn; i++)
                        ddteta[i] = 6 * a3 + qy * x[i] / EI;

                }
                else if (ElType == "Beam2t_hl")
                {   //Timoshenko beam element - left hinge

                    double v1, t1, v2, t2;
                    double a0, a1, a2, a3;
                    double ph, phs, phss;

                    ph = 12 * EI / (GAK * Math.Pow(L, 2));
                    phs = 1 / (1 + ph);
                    phss = EI / GAK;

                    v1 = m[1];
                    t2 = m[5];
                    v2 = m[4];
                    t1 = -(1 + ph) / (4 * (1 + ph / 4)) * (6 * phs / L * v1 + phs * 2 * (1 - ph / 2) * t2 - 6 * phs / L * v2);

                    a0 = v1;
                    a1 = phs * (-1 / L * ph * v1 + (1 + ph / 2) * t1 + 1 / L * ph * v2 - ph / 2 * t2);
                    a2 = 0;
                    a3 = phs * (2 * v1 / Math.Pow(L, 3) + t1 / Math.Pow(L, 2) - 2 * v2 / Math.Pow(L, 3) + t2 / Math.Pow(L, 2));


                    for (int i = 0; i < nn; i++)
                        v[i] = Math.Pow(x[i], 3) * a3 + Math.Pow(x[i], 2) * a2 + x[i] * a1 + a0 + qy / (24 * EI) * Math.Pow(x[i], 4) - qy / (2 * GAK) * Math.Pow(x[i], 2);

                    for (int i = 0; i < nn; i++)
                        dv[i] = 3 * Math.Pow(x[i], 2) * a3 + 2 * x[i] * a2 + a1 + qy * Math.Pow(x[i], 3) / (6 * EI) - qy * x[i] / GAK;

                    for (int i = 0; i < nn; i++)
                        teta[i] = 3 * (Math.Pow(x[i], 2) + 2 * alfa) * a3 + a1 + qy * Math.Pow(x[i], 3) / (6 * EI);

                    for (int i = 0; i < nn; i++)
                        dteta[i] = 6 * x[i] * a3 + qy * Math.Pow(x[i], 2) / (2 * EI);

                    for (int i = 0; i < nn; i++)
                        ddteta[i] = 6 * a3 + qy * x[i] / EI;

                }
                else if (ElType == "Beam2t_hr")
                {   //Timoshenko beam element - right hinge

                    double v1, t1, v2, t2;
                    double a0, a1, a2, a3;
                    double ph, phs, phss;

                    ph = 12 * EI / (GAK * Math.Pow(L, 2));
                    phs = 1 / (1 + ph);
                    phss = EI / GAK;



                    for (int i = nn - 1; i >= 0; i--)
                        x[nn - 1 - i] = i * ll;



                    v1 = m[4];
                    t2 = -m[2];
                    v2 = m[1];
                    t1 = -(1 + ph) / (4 * (1 + ph / 4)) * (6 * phs / L * v1 + phs * 2 * (1 - ph / 2) * t2 - 6 * phs / L * v2);

                    a0 = v1;
                    a1 = phs * (-1 / L * ph * v1 + (1 + ph / 2) * t1 + 1 / L * ph * v2 - ph / 2 * t2);
                    a2 = 0;
                    a3 = phs * (2 * v1 / Math.Pow(L, 3) + t1 / Math.Pow(L, 2) - 2 * v2 / Math.Pow(L, 3) + t2 / Math.Pow(L, 2));



                    for (int i = 0; i < nn; i++)
                        v[i] = Math.Pow(x[i], 3) * a3 + Math.Pow(x[i], 2) * a2 + x[i] * a1 + a0 + qy / (24 * EI) * Math.Pow(x[i], 4) - qy / (2 * GAK) * Math.Pow(x[i], 2);

                    for (int i = 0; i < nn; i++)
                        dv[i] = 3 * Math.Pow(x[i], 2) * a3 + 2 * x[i] * a2 + a1 + qy * Math.Pow(x[i], 3) / (6 * EI) - qy * x[i] / GAK;

                    for (int i = 0; i < nn; i++)
                        teta[i] = 3 * (Math.Pow(x[i], 2) + 2 * alfa) * a3 + a1 + qy * Math.Pow(x[i], 3) / (6 * EI);

                    for (int i = 0; i < nn; i++)
                        dteta[i] = 6 * x[i] * a3 + qy * Math.Pow(x[i], 2) / (2 * EI);

                    for (int i = 0; i < nn; i++)
                        ddteta[i] = 6 * a3 + qy * x[i] / EI;

                }


                for (int i = 0; i < nn; i++)
                    lc.Du[el].Add(u[i]);

                for (int i = 0; i < nn; i++)
                    lc.Dv[el].Add(v[i]);

                for (int i = 0; i < nn; i++)
                    lc.N[el].Add(EA * du[i]);

                for (int i = 0; i < nn; i++)
                    lc.M[el].Add(EI * dteta[i]);

                for (int i = 0; i < nn; i++)
                    lc.V[el].Add(EI * ddteta[i]);


            }

            /*
            //Add nodal load (nn, fx, fy, m)
            public void addnLoad(int nnode, nodeLoad nl)
            {

                Nodes[nnode].nl.Add(nl);

            }
            */

            //Check if the node is new
            private int isNew(Node2d nod)
            {

                int num = 0;

                for (int i = 0; i < Nodes.Count; i++)
                    if ((Nodes[i].x == nod.x) && (Nodes[i].y == nod.y))
                        num = Nodes[i].num;

                return num;
            }

            //Add node
            public int addNode(Node2d nod)
            {

                int nn = isNew(nod);

                if (nn == 0)
                {
                    nod.num = Nodes.Count + 1;
                    Nodes.Add(nod);
                }
                else
                    nod.num = nn;

                return nod.num;
            }

            //Add element
            public void addEl2d(El2d e2d)
            {

                int n1, n2;

                Elements.Add(e2d);

                n1 = addNode(e2d.Nodes[0]);
                n2 = addNode(e2d.Nodes[1]);

                e2d.Nodes[0].num = n1;
                e2d.Nodes[1].num = n2;


            }

            //Assemble global stiffness matrix
            public void assem(double[][] Ke, double[] fe, List<int> edofe)
            {

                int n = edofe.Count;

                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++)
                        K[edofe[i] - 1][edofe[j] - 1] = K[edofe[i] - 1][edofe[j] - 1] + Ke[i][j];

                for (int i = 0; i < n; i++)
                    f[edofe[i] - 1] = f[edofe[i] - 1] + fe[i];

            }

            //Solve system with LU factorization
            public void Solve()
            {


                LUdecomp(); //LU decomposition


                int i, j;

                int maxOrder = l.Length;

                disp = new double[maxOrder];

                f0 = Matrix.ArrayDuplicate(f);


                for (j = 1; j < maxOrder; j++)
                {
                    for (i = 0; i < j; i++)
                    {
                        f0[j] = f0[j] - f0[i] * l[j][i];
                    }
                }

                disp[maxOrder - 1] = f0[maxOrder - 1] / u[maxOrder - 1][maxOrder - 1];

                for (j = maxOrder - 2; j >= 0; j--)
                {
                    disp[j] = f0[j];
                    for (i = maxOrder - 1; i > j; i--)
                    {
                        disp[j] = disp[j] - disp[i] * u[j][i];
                    }
                    disp[j] = disp[j] / u[j][j];
                }
            }

            //Set boundary condition (BCs-Penalty approach)
            public void setBC(LoadCase lc)
            {

                int i, j, k;
                List<int> edofe;

                ndof = Nodes[Nodes.Count - 1].ndof3;

                K = new double[ndof][];
                for (i = 0; i < ndof; i++)
                    K[i] = new double[ndof];

                f = new double[ndof];


                //Initialize internal forces
                lc.Init(Elements.Count);

                double[][] fle = lc.genFle(Elements);

                //Assembling matrix before setting the boundary conditions
                for (i = 0; i < Elements.Count; i++)
                {

                    edofe = new List<int> { Elements[i].Nodes[0].ndof1, Elements[i].Nodes[0].ndof2, Elements[i].Nodes[0].ndof3, Elements[i].Nodes[1].ndof1, Elements[i].Nodes[1].ndof2, Elements[i].Nodes[1].ndof3 };

                    edof.Add(edofe);

                    assem(Elements[i].K, fle[i], edofe); //lc.fle

                }


                int nk = K.Length;

                double Kmax = 0;

                Kbc = Matrix.MatrixDuplicate(K);

                for (i = 0; i < nk; i++)
                    for (j = i; j < nk; j++)
                        if (Math.Abs(Kbc[i][j]) > Kmax)
                            Kmax = Math.Abs(Kbc[i][j]);

                double C = 1000000 * Kmax; //penalty approach -> 10^6


                for (i = 0; i < Nodes.Count; i++)
                    if (Nodes[i].bc != null)
                        for (j = 0; j < nk; j++)
                            if (((Nodes[i].bc.vx == 1) && (Nodes[i].ndof1 - 1 == j)) || ((Nodes[i].bc.vy == 1) && (Nodes[i].ndof2 - 1 == j)) || ((Nodes[i].bc.rr == 1) && (Nodes[i].ndof3 - 1 == j)))
                                Kbc[j][j] = Kbc[j][j] + C;


                for (i = 0; i < lc.nl.Count; i++)
                {
                    f[Nodes[lc.nl[i].nn].ndof1 - 1] = f[Nodes[lc.nl[i].nn].ndof1 - 1] + lc.nl[i].fx;
                    f[Nodes[lc.nl[i].nn].ndof2 - 1] = f[Nodes[lc.nl[i].nn].ndof2 - 1] + lc.nl[i].fy;
                    f[Nodes[lc.nl[i].nn].ndof3 - 1] = f[Nodes[lc.nl[i].nn].ndof3 - 1] + lc.nl[i].m;
                }

            }

            //Calculate reactions
            public void calcReact()
            {

                int i, j, k;
                int nk = K.Length;

                React = new double[nk];

                for (i = 0; i < nk; i++)
                {

                    React[i] = 0;
                    for (k = 0; k < nk; k++)
                        React[i] = React[i] + K[i][k] * disp[k];

                    React[i] = React[i] - f[i];

                }

            }

            //LU decomposition
            public void LUdecomp()
            {

                int i, j, k;
                int maxOrder = Kbc.Length;
                l = new double[maxOrder][];
                u = new double[maxOrder][];

                for (i = 0; i < maxOrder; i++)
                {

                    l[i] = new double[maxOrder];
                    u[i] = new double[maxOrder];

                }

                for (i = 0; i < maxOrder; i++)
                    l[i][i] = 1.0;

                for (i = 0; i < maxOrder; i++)
                {
                    if (Math.Abs(Kbc[i][i]) < 1E-10)
                        Matrix.SwitchRows(i, f, Kbc);
                }
                for (j = 0; j < maxOrder; j++)
                {
                    for (i = 0; i < maxOrder; i++)
                    {
                        if (i >= j)
                        {
                            u[j][i] = Kbc[j][i];
                            for (k = 0; k < j; k++)
                                u[j][i] = u[j][i] - u[k][i] * l[j][k];
                        }
                        if (i > j)
                        {
                            l[i][j] = Kbc[i][j];
                            for (k = 0; k < j; k++)
                                l[i][j] = l[i][j] - u[k][j] * l[i][k];
                            l[i][j] = l[i][j] / u[j][j];
                        }
                    }
                }
            }

        }

        static void Main(string[] args)
        {

            //Create structure
            Str2d str = new Str2d();

            //Create nodes
            str.addNode(new Node2d(1, 0, 0));
            str.addNode(new Node2d(2, 0, 3000));
            str.addNode(new Node2d(3, 3000, 3000));
            str.addNode(new Node2d(4, 3000, 0));

            //Shear deformation
            bool shearDef = false;
            double G = 81000; //Steel

            if (!shearDef)
                G *= 1000000;

            //Define material and cross section
            Material mat = new Material(210000, G); //G=81000
            Section sec = new Section(10000, 8333333, 0.8);

            //Create elements
            str.addEl2d(new El2d(0, "Beam2t", str.Nodes[0], str.Nodes[1], mat, sec));
            str.addEl2d(new El2d(1, "Beam2t_hr", str.Nodes[1], str.Nodes[2], mat, sec));
            str.addEl2d(new El2d(2, "Beam2t", str.Nodes[2], str.Nodes[3], mat, sec));

            //Define boundary conditions
            str.Nodes[0].bc = new BC(0, 1, 1, 1);
            str.Nodes[3].bc = new BC(3, 1, 1, 1);



            LoadCase lc = new LoadCase("Perm");

            //Create node load
            //nodeLoad nlp = new nodeLoad(1, 10000, 0, 0); //N
            //lc.addNodeLoad(nlp);

            //Create uniform load
            unifLoad ulp = new unifLoad(1, 0, -1); //N/mm -> kN/m
            lc.addUnifLoad(ulp);




            //Set boundary conditions
            str.setBC(lc);

            //Solve system
            str.Solve();

            //Calculate reactions
            str.calcReact();

            //Extract end-displacements for each element
            str.extract();

            //Calculate internal forces
            str.IntForces(lc, 0, 2);
            str.IntForces(lc, 1, 11);
            str.IntForces(lc, 2, 2);

        }

    }
}
