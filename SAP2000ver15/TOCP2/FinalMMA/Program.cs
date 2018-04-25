using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FinalMMA
{
    class Program
    {
        //METHODS
        public void ReadModel(ref SAP2000v15.cSapModel Model, out string[] PointName, out int PointNumb, out string[] SolidName, out int SolidNumb, out string[] FrameName, out int FrameNumb, out string[] AreaName, out int AreaNumb, out string[] LoadName, out int LoadNumb)
        {
            int i;
            int ret;
            //Read points and change their local axes 
            PointName = new string[0];
            PointNumb = 0;
            double a = 0.0;
            double b = 0.0;
            double c = 0.0;
            ret = Model.PointObj.GetNameList(ref PointNumb, ref PointName);
            Console.WriteLine("{0} points", PointNumb);
            for (i = 0; i <= PointNumb - 1; i++)
            {
                ret = Model.PointObj.SetLocalAxes(PointName[i], a, b, c, SAP2000v15.eItemType.Object);
            }
            //Read solid objects and change their local axes 
            SolidName = new string[0];
            SolidNumb = 0;
            ret = Model.SolidObj.GetNameList(ref SolidNumb, ref SolidName);
            Console.WriteLine("{0} solid objects", SolidNumb);
            for (i = 0; i <= SolidNumb - 1; i++)
            {
                ret = Model.SolidObj.SetLocalAxes(SolidName[i], a, b, c, SAP2000v15.eItemType.Object);
            }
            //Read frame objects and change their local axes 
            FrameName = new string[0];
            FrameNumb = 0;
            ret = Model.FrameObj.GetNameList(ref FrameNumb, ref FrameName);
            Console.WriteLine("{0} Frame objects ", FrameNumb);
            for (i = 0; i <= FrameNumb - 1; i++)
            {
                ret = Model.FrameObj.SetLocalAxes(FrameName[i], 0.0, SAP2000v15.eItemType.Object);
            }
            //Read areaobject and change their local axes 
            AreaName = new string[0];
            AreaNumb = 0;
            ret = Model.AreaObj.GetNameList(ref AreaNumb, ref AreaName);
            Console.WriteLine("{0} Area objects ", AreaNumb);
            for (i = 0; i <= AreaNumb - 1; i++)
            {
                ret = Model.AreaObj.SetLocalAxes(AreaName[i], 0.0, SAP2000v15.eItemType.Object);
            }
            //Read loadcases
            LoadNumb = 0;
            LoadName = new string[0];
            ret = Model.LoadCases.GetNameList(ref LoadNumb, ref LoadName);
        }
        public void CreateMatSolid(ref SAP2000v15.cSapModel Model, string Region, string Standard, string Grade)
        {
            int i;
            int ret;
            string temp_name;
            string[] MatName;
            MatName = new string[101];
            string temp_mat = "";
            double e = 0.0;
            double amat = 0.0;
            double u = 0.0;
            double g = 0.0;
            //Create 101 materials with names OptiMat0, OptiMat1..
            for (i = 0; i <= 100; i++)
            {
                temp_name = "OptiMat" + Convert.ToString(i);
                ret = Model.PropMaterial.AddMaterial(ref temp_mat, SAP2000v15.eMatType.MATERIAL_CONCRETE, Region, Standard, Grade, temp_name);
                MatName[i] = temp_mat;
            }
            //Change of E in each material
            ret = Model.PropMaterial.GetMPIsotropic(MatName[0], ref e, ref u, ref amat, ref g);
            ret = Model.PropMaterial.SetMPIsotropic(MatName[0], e * 0.0000001, u, amat);
            double metr = 0.01;
            for (i = 1; i <= 100; i++)
            {
                ret = Model.PropMaterial.SetMPIsotropic(MatName[i], e * metr, u, amat);
                metr = metr + 0.01;
                metr = Math.Round(metr, 2);
            }
            //Create 101 solid properties with names OptiSolidProp0.00, OptiSolidProp0.01..
            metr = 0.00;
            for (i = 0; i <= 100; i++)
            {
                temp_name = "OptiSolidProp" + string.Format("{0:N2}", metr);
                ret = Model.PropSolid.SetProp(temp_name, MatName[i], 0.0, 0.0, 0.0, true, -1, "", "");
                metr = Math.Round(metr, 2) + 0.01;
                metr = Math.Round(metr, 2);
            }
        }   
        public void CreateOptiSolidName(ref SAP2000v15.cSapModel Model, ref bool testbool, int arstoix, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double diakr, out int OptiSolidNumb, out string[] OptiSolidName)
        {
            int ret;
            double x_double = 0.0;
            double y_double = 0.0;
            double z_double = 0.0;
            OptiSolidName = new string[arstoix];
            int count = 0;
            int NumbSel = 0;
            int[] ObjTypeSel;
            ObjTypeSel = new int[1];
            string[] ObjNameSel;
            ObjNameSel = new string[1];
            string voith = "";
            for (x_double =  Math.Round(xmax,3); x_double >= Math.Round(xmin,3) + Math.Round(diakr/2, 3); x_double = Math.Round(x_double, 3) - Math.Round(diakr, 3))
            {
                for (z_double =  Math.Round(zmax,3); z_double >= Math.Round(zmin,3) + Math.Round(diakr/2, 3); z_double = Math.Round(z_double, 3) - Math.Round(diakr, 3))
                {
                    for (y_double =  Math.Round(ymax,3); y_double >= Math.Round(ymin,3) + Math.Round(diakr/2, 3); y_double = Math.Round(y_double, 3) - Math.Round(diakr, 3))
                    {
                        ret = Model.SelectObj.ClearSelection();
                        ret = Model.SelectObj.CoordinateRange(Math.Round(x_double,3) - Math.Round(diakr,3), Math.Round(x_double,3), Math.Round(y_double,3) - Math.Round(diakr,3), Math.Round(y_double,3), Math.Round(z_double,3) - Math.Round(diakr,3), Math.Round(z_double,3), false, "Global", false, false, false, false, true, false);
                        if (ret != 0)
                        {
                            Console.WriteLine("Error: while creating OptiSolidName.");
                        }
                        ret = Model.SelectObj.GetSelected(ref NumbSel, ref ObjTypeSel, ref ObjNameSel);
                        if (ret != 0)
                        {
                            Console.WriteLine("Error: while creating OptiSolidName.");
                        }
                        if (ObjNameSel[0] == voith)
                        {
                            Console.WriteLine("Error: while creating OptiSolidName.");
                        }
                        count = count + 1;
                        OptiSolidName[count - 1] = ObjNameSel[0];
                        voith = ObjNameSel[0];
                    }
                }
            }
            OptiSolidNumb = count;
            if (OptiSolidNumb != arstoix)
            {
                Console.WriteLine("Error: while creating OptiSolidName. OptiSolidNumb is {0} while arstoix is {1}. check if the program didnt read all the elements", OptiSolidNumb, arstoix);
                testbool = false;
            }
        }
        public void CreateMapDist(ref SAP2000v15.cSapModel Model, int OptiSolidNumb, int nlx, int nly, int nlz, double diakr, double rmin, out int[,] Map, out double[,] Dist)
        {
            int i;
            int j;
            int x_int;
            int y_int;
            int z_int;
            double x_double_max;
            double y_double_max;
            double z_double_max;
            double x_double_min;
            double y_double_min;
            double z_double_min;
            double x_double;
            double y_double;
            double z_double;
            Map = new int[OptiSolidNumb, 27];
            Dist = new double[OptiSolidNumb, 27];
            int e_int;
            double e_double;
            double H;
            for (i = 0; i <= OptiSolidNumb - 1; i++)
            {
                for (j = 0; j <= 26; j++)
                {
                    Map[i, j] = 0;
                    Dist[i, j] = 0.0;
                }
            }
            for (x_int = 1; x_int <= nlx; x_int++)
            {
                for (z_int = 1; z_int <= nlz; z_int++)
                {
                    for (y_int = 1; y_int <= nly; y_int++)
                    {
                        e_int = (x_int - 1) * nlz * nly + (z_int - 1) * nly + y_int;
                        x_double_min = Math.Truncate(Math.Max(x_int - rmin, 0.0)) + 1;
                        z_double_min = Math.Truncate(Math.Max(z_int - rmin, 0.0)) + 1;
                        y_double_min = Math.Truncate(Math.Max(y_int - rmin, 0.0)) + 1;
                        x_double_max = Math.Min(Math.Truncate(x_int + rmin), nlx);
                        z_double_max = Math.Min(Math.Truncate(z_int + rmin), nlz);
                        y_double_max = Math.Min(Math.Truncate(y_int + rmin), nly);
                        j = 0;
                        for (x_double = x_double_min; x_double <= x_double_max; x_double++)
                        {
                            for (z_double = z_double_min; z_double <= z_double_max; z_double++)
                            {
                                for (y_double = y_double_min; y_double <= y_double_max; y_double++)
                                {
                                    e_double = (x_double - 1) * nlz * nly + (z_double - 1) * nly + y_double;
                                    H = (rmin - Math.Sqrt(Math.Pow(x_int - Convert.ToInt32(x_double), 2) + Math.Pow(y_int - Convert.ToInt32(y_double), 2) + Math.Pow(z_int - Convert.ToInt32(z_double), 2))) * diakr;
                                    Map[e_int - 1, j] = Convert.ToInt32(e_double);
                                    Dist[e_int - 1, j] = Math.Max(H, 0);
                                    j = j + 1;
                                }
                            }
                        }
                    }
                }
            }
        }
        public void UpdateModel(ref SAP2000v15.cSapModel Model, string[] OptiSolidName, int OptiSolidNumb, double[] xkfil)
        {
            int i;
            int ret;
            string temp_string;
            double temp_double;
            //update model based on array xkfil           
            ret = Model.SetModelIsLocked(false);
            for (i = 0; i <= OptiSolidNumb - 1; i++)
            {
                temp_double = Math.Round(xkfil[i] * xkfil[i] * xkfil[i], 2);
                temp_string = "OptiSolidProp" + string.Format("{0:N2}", temp_double);
                ret = Model.SolidObj.SetProperty(OptiSolidName[i], temp_string, 0);
            }
        }
        public double CalcCompliance(ref SAP2000v15.cSapModel Model, ref bool testbool, string[] SolidName, int SolidNumb, string[] FrameName, int FrameNumb, string[] AreaName, int AreaNumb)
        {
            //Declaration for Joint Displacements
            int nresult = 0;
            string[] LoadCase;
            LoadCase = new string[1];
            string[] StepType;
            StepType = new string[1];
            double[] StepNum;
            StepNum = new double[1];
            double[] U1;
            U1 = new double[1];
            double[] U2;
            U2 = new double[1];
            double[] U3;
            U3 = new double[1];
            double[] R1;
            R1 = new double[1];
            double[] R2;
            R2 = new double[1];
            double[] R3;
            R3 = new double[1];
            string[] Obje;
            Obje = new string[1];
            string[] elm;
            elm = new string[1];

            //Declaration for solid joint force
            int nresultS = 0;
            string[] ObjeS;
            ObjeS = new string[8];
            string[] elmS;
            elmS = new string[8];
            string[] LoadCaseS;
            LoadCaseS = new string[8];
            string[] StepTypeS;
            StepTypeS = new string[8];
            double[] StepNumS;
            StepNumS = new double[8];
            double[] F1S;
            F1S = new double[8];
            double[] F2S;
            F2S = new double[8];
            double[] F3S;
            F3S = new double[8];
            double[] M1S;
            M1S = new double[8];
            double[] M2S;
            M2S = new double[8];
            double[] M3S;
            M3S = new double[8];
            string[] pointelmS;
            pointelmS = new string[8];

            //Declaration for frame joint forces
            int nresultF = 0;
            string[] LoadCaseF;
            LoadCaseF = new string[2];
            string[] StepTypeF;
            StepTypeF = new string[2];
            double[] StepNumF;
            StepNumF = new double[2];
            string[] pointelmF;
            pointelmF = new string[2];
            double[] F1F;
            F1F = new double[2];
            double[] F2F;
            F2F = new double[2];
            double[] F3F;
            F3F = new double[2];
            double[] M1F;
            M1F = new double[2];
            double[] M2F;
            M2F = new double[2];
            double[] M3F;
            M3F = new double[2];
            string[] ObjeF;
            ObjeF = new string[2];
            string[] elmF;
            elmF = new string[2];

            //Declaration for AreaJointForceShell and AreaJointForcePlane
            int nresultA = 0;
            string[] LoadCaseA;
            LoadCaseA = new string[4];
            string[] StepTypeA;
            StepTypeA = new string[4];
            double[] StepNumA;
            StepNumA = new double[4];
            string[] pointelmA;
            pointelmA = new string[4];
            double[] F1A;
            F1A = new double[4];
            double[] F2A;
            F2A = new double[4];
            double[] F3A;
            F3A = new double[4];
            double[] M1A;
            M1A = new double[4];
            double[] M2A;
            M2A = new double[4];
            double[] M3A;
            M3A = new double[4];
            string[] ObjeA;
            ObjeA = new string[4];
            string[] elmA;
            elmA = new string[4];
            int proptype = 0;
            string propname = "";

            //COMPLIANCE
            int i;
            int j;
            int ret;
            double comp = 0.0;
            if (SolidNumb != 0)
            {
                for (i = 0; i <= SolidNumb - 1; i++)
                {
                    ret = Model.Results.SolidJointForce(SolidName[i], SAP2000v15.eItemTypeElm.ObjectElm, ref nresultS, ref ObjeS, ref elmS, ref pointelmS, ref LoadCaseS, ref StepTypeS, ref StepNumS, ref F1S, ref F2S, ref F3S, ref M1S, ref M2S, ref M3S);
                    if (nresultS == 8)
                    {
                        for (j = 0; j <= nresultS - 1; j++)
                        {
                            ret = Model.Results.JointDispl(pointelmS[j], SAP2000v15.eItemTypeElm.ObjectElm, ref nresult, ref Obje, ref elm, ref LoadCase, ref StepType, ref StepNum, ref U1, ref U2, ref U3, ref R1, ref R2, ref R3);
                            comp = comp + U1[0] * F1S[j] + U2[0] * F2S[j] + U3[0] * F3S[j] + R1[0] * M1S[j] + R2[0] * M2S[j] + R3[0] * M3S[j];
                        }
                    }
                    else
                    {
                        Console.WriteLine("Error: while calculating the compliance of a solid element (Method CalcCompliance). check nresultS");
                        Console.WriteLine("number of results: {0}", nresultS);
                        testbool = false;
                    }
                }
            }
            if (FrameNumb != 0)
            {
                for (i = 0; i <= FrameNumb - 1; i++)
                {
                    ret = Model.Results.FrameJointForce(FrameName[i], SAP2000v15.eItemTypeElm.ObjectElm, ref nresultF, ref ObjeF, ref elmF, ref pointelmF, ref LoadCaseF, ref StepTypeF, ref StepNumF, ref F1F, ref F2F, ref F3F, ref M1F, ref M2F, ref M3F);
                    if (nresultF == 2)
                    {
                        for (j = 0; j <= nresultF - 1; j++)
                        {
                            ret = Model.Results.JointDispl(pointelmF[j], SAP2000v15.eItemTypeElm.ObjectElm, ref nresult, ref Obje, ref elm, ref LoadCase, ref StepType, ref StepNum, ref U1, ref U2, ref U3, ref R1, ref R2, ref R3);
                            comp = comp + U1[0] * F1F[j] + U2[0] * F2F[j] + U3[0] * F3F[j] + R1[0] * M1F[j] + R2[0] * M2F[j] + R3[0] * M3F[j];
                        }
                    }
                    else
                    {
                        Console.WriteLine("Error: while calculating the compliance of a Frame element (Method CalcCompliance). check nresultF");
                        Console.WriteLine("number of results: {0}", nresultF);
                        testbool = false;
                    }
                }
            }
            if (AreaNumb != 0)
            {
                for (i = 0; i <= AreaNumb - 1; i++)
                {
                    ret = Model.AreaObj.GetProperty(AreaName[i], ref propname);
                    ret = Model.PropArea.GetType(propname, ref proptype); //1 = Shell 2 = Plane 3 = Asolid               
                    if (proptype == 1)
                    {
                        ret = Model.Results.AreaJointForceShell(AreaName[i], SAP2000v15.eItemTypeElm.ObjectElm, ref nresultA, ref ObjeA, ref elmA, ref pointelmA, ref LoadCaseA, ref StepTypeA, ref StepNumA, ref F1A, ref F2A, ref F3A, ref M1A, ref M2A, ref M3A);
                        if (nresultA == 4)
                        {
                            for (j = 0; j <= nresultA - 1; j++)
                            {
                                ret = Model.Results.JointDispl(pointelmA[j], SAP2000v15.eItemTypeElm.ObjectElm, ref nresult, ref Obje, ref elm, ref LoadCase, ref StepType, ref StepNum, ref U1, ref U2, ref U3, ref R1, ref R2, ref R3);

                                comp = comp + U1[0] * F1A[j] + U2[0] * F2A[j] + U3[0] * F3A[j] + R1[0] * M1A[j] + R2[0] * M2A[j] + R3[0] * M3A[j];
                            }
                        }
                        else
                        {
                            Console.WriteLine("Error: while calculating the compliance of an area-shell element (Method CalcCompliance). check nresultA");
                            Console.WriteLine("number of results: {0}", nresultA);
                            testbool = false;
                        }
                    }
                    else if (proptype == 2)
                    {
                        ret = Model.Results.AreaJointForcePlane(AreaName[i], SAP2000v15.eItemTypeElm.ObjectElm, ref nresultA, ref ObjeA, ref elmA, ref pointelmA, ref LoadCaseA, ref StepTypeA, ref StepNumA, ref F1A, ref F2A, ref F3A, ref M1A, ref M2A, ref M3A);
                        if (nresultA == 4)
                        {
                            for (j = 0; j <= nresultA - 1; j++)
                            {
                                ret = Model.Results.JointDispl(pointelmA[j], SAP2000v15.eItemTypeElm.ObjectElm, ref nresult, ref Obje, ref elm, ref LoadCase, ref StepType, ref StepNum, ref U1, ref U2, ref U3, ref R1, ref R2, ref R3);
                                comp = comp + U1[0] * F1A[j] + U2[0] * F2A[j] + U3[0] * F3A[j] + R1[0] * M1A[j] + R2[0] * M2A[j] + R3[0] * M3A[j];
                            }
                        }
                        else
                        {
                            Console.WriteLine("Error: while calculating the compliance of an area-plane element(Method CalcCompliance). check nresultA");
                            Console.WriteLine("number of results: {0}", nresultA);
                            testbool = false;
                        }
                    }
                    else
                    {
                        Console.WriteLine("Error: while calculating the compliance of an area element(Method CalcCompliance). check proptype ");
                        Console.WriteLine("proptype = {0}", proptype);
                        testbool = false;
                    }
                }
            }
            return comp;
        }
        public void CalcDerivative(ref SAP2000v15.cSapModel Model, ref bool testbool, string[] OptiSolidName, int OptiSolidNumb, double[] xkfil, out double[] der, out double comp2)
        {
            //declaration for solid joint force
            int nresultS = 0;
            string[] ObjeS;
            ObjeS = new string[8];
            string[] elmS;
            elmS = new string[8];
            string[] LoadCaseS;
            LoadCaseS = new string[8];
            string[] StepTypeS;
            StepTypeS = new string[8];
            double[] StepNumS;
            StepNumS = new double[8];
            double[] F1S;
            F1S = new double[8];
            double[] F2S;
            F2S = new double[8];
            double[] F3S;
            F3S = new double[8];
            double[] M1S;
            M1S = new double[8];
            double[] M2S;
            M2S = new double[8];
            double[] M3S;
            M3S = new double[8];
            string[] pointelmS;
            pointelmS = new string[8];

            //declaration for Joint Displacements
            int nresult = 0;
            string[] LoadCase;
            LoadCase = new string[1];
            string[] StepType;
            StepType = new string[1];
            double[] StepNum;
            StepNum = new double[1];
            double[] U1;
            U1 = new double[1];
            double[] U2;
            U2 = new double[1];
            double[] U3;
            U3 = new double[1];
            double[] R1;
            R1 = new double[1];
            double[] R2;
            R2 = new double[1];
            double[] R3;
            R3 = new double[1];
            string[] Obje;
            Obje = new string[1];
            string[] elm;
            elm = new string[1];

            //rest of declarations
            int i;
            int j;
            int ret;
            der = new double[OptiSolidNumb];
            double compder;
            comp2 = 0.0;
            for (i = 0; i <= OptiSolidNumb - 1; i++)
            {
                compder = 0.0;
                ret = Model.Results.SolidJointForce(OptiSolidName[i], SAP2000v15.eItemTypeElm.ObjectElm, ref nresultS, ref ObjeS, ref elmS, ref pointelmS, ref LoadCaseS, ref StepTypeS, ref StepNumS, ref F1S, ref F2S, ref F3S, ref M1S, ref M2S, ref M3S);
                if (nresultS == 8)
                {
                    for (j = 0; j <= nresultS - 1; j++)
                    {
                        ret = Model.Results.JointDispl(pointelmS[j], SAP2000v15.eItemTypeElm.ObjectElm, ref nresult, ref Obje, ref elm, ref LoadCase, ref StepType, ref StepNum, ref U1, ref U2, ref U3, ref R1, ref R2, ref R3);
                                                                                             
                        compder = compder + (U1[0] * F1S[j]) + (U2[0] * F2S[j]) + (U3[0] * F3S[j]) + (R1[0] * M1S[j]) + (R2[0] * M2S[j]) + (R3[0] * M3S[j]);

                    }
                }
                else
                {
                    Console.WriteLine("Error: while calculating the compliance of a solid element (Method CalcDerivatives). check nresultS");
                    Console.WriteLine("number of results: {0}", nresultS);
                    testbool = false;
                }
                der[i] = -3 * compder / Math.Max(xkfil[i], Math.Pow(10,-3));
                if (der[i] > 0)
                {
                    Console.WriteLine("Error: positive derivative = {0} in solid element {1}", der[i], OptiSolidName[i]);
                    testbool = false;
                }
                comp2 = comp2 + compder;
            }
        }
        public double[] FilterDer(ref bool testbool, int OptiSolidNumb, int[,] Map, double[,] Dist, double[] xkfil, double[] der)
        {
            int i;
            int j;
            double SHXC;
            double SH;
            double[] derfil;
            derfil = new double[OptiSolidNumb];
            for (i = 0; i <= OptiSolidNumb - 1; i++)
            {
                SHXC = 0.0;
                SH = 0.0;
                for (j = 0; j <= 26; j++)
                {
                    if ((Map[i, j] - 1) >= 0)
                    {
                        SHXC = SHXC + Dist[i, j] * xkfil[Map[i, j] - 1] * der[Map[i, j] - 1];
                        SH = SH + Dist[i, j];
                    }
                    else if ((Map[i, j] - 1) != -1)
                    {
                        Console.WriteLine("Error:  Map[i,j]-1 is {0}", Map[i, j] - 1); 
                        testbool = false;
                    }
                }
                derfil[i] = SHXC / (SH * Math.Max(xkfil[i],Math.Pow(10,-3)));
            }
            return derfil;
        }
        public void OptimalityCriteria(int OptiSolidNumb, double[] xk, double[] derfil, int[,] Map, double[,] Dist, double V, out double[] xnew, out double[] xnewfil)
        {
            int i;
            int j;
            double SH;
            double lamda;
            double SHX;
            double metrx;
            double al;
            double bl;
            double[] be;
            be = new double[OptiSolidNumb];
            xnew = new double[OptiSolidNumb];
            xnewfil = new double[OptiSolidNumb];
            bool testlamda = false;
            for (i = 0; i <= OptiSolidNumb - 1; i++)
            {
                be[i] = -xk[i] * xk[i] * derfil[i];
            }
           testlamda = false;
           al = 0.0;
           bl = Math.Pow(10, 10);                    
           do
            {
                lamda = (al + bl) / 2;
                for (i = 0; i <= OptiSolidNumb - 1; i++)
                {
                    xnew[i] = Math.Max(0, Math.Max(xk[i] - 0.2, Math.Min(1, Math.Min(xk[i] + 0.2, Math.Sqrt(be[i] / lamda)))));
                }
                metrx = 0.0;
                for (i = 0; i <= OptiSolidNumb - 1; i++)
                {
                    SHX = 0.0;
                    SH = 0.0;
                    for (j = 0; j <= 26; j++)
                    {
                        if ((Map[i, j] - 1) != -1)
                        {
                            SHX = SHX + Dist[i, j] * xnew[Map[i, j] - 1];
                            SH = SH + Dist[i, j];
                        }
                    }
                    xnewfil[i] = SHX / SH;
                    metrx = metrx + xnewfil[i];
                }
                if ((bl - al) / (al + bl) < Math.Pow(10, -3))
                {
                    testlamda = true;
                }
                else if (metrx > V)
                {
                    al = lamda;
                }
                else
                {
                    bl = lamda;
                }
            } while (testlamda == false);
        }     
        public void MMA(ref bool testbool, int epanal, int OptiSolidNumb, double[] xk, double[] xkm1, double[] xkm2, double[] derfil, int[,] Map, double[,] Dist, double V, out double[] xnew, out double[] xnewfil)
        {
            int i;
            int j;
            double al;
            double bl;
            bool testlamda;
            double gamma;
            double lamda;
            double SHX;
            double SH;
            double metrx;
            xnew = new double[OptiSolidNumb];
            xnewfil = new double[OptiSolidNumb];
            //declaration L U
            double[] Lk;
            Lk = new double[OptiSolidNumb];
            double[] Ukm1;
            Ukm1 = new double[OptiSolidNumb];
            double[] Uk;
            Uk = new double[OptiSolidNumb];
            //declaration a b
            double[] ak;
            ak = new double[OptiSolidNumb];
            double[] bk;
            bk = new double[OptiSolidNumb];
            //declaration pk qk
            double[] pk;
            pk = new double[OptiSolidNumb];
            double[] qk;
            qk = new double[OptiSolidNumb];
            if (epanal >= 3)
            {
                for (i = 0; i <= OptiSolidNumb - 1; i++)
                {
                    if ((xk[i] - xkm1[i]) * (xkm1[i] - xkm2[i]) < 0)
                    {
                        gamma = 0.7;
                    }
                    else if ((xk[i] - xkm1[i]) * (xkm1[i] - xkm2[i]) > 0)
                    {
                        gamma = 1.2;
                    }
                    else
                    {
                        gamma = 1.0;
                    }
                    Lk[i] = xk[i] - 0.5 * gamma;
                    Uk[i] = xk[i] + 0.5 * gamma;

                    ak[i] = 0.9 * Lk[i] + 0.1 * xk[i];
                    bk[i] = 0.9 * Uk[i] + 0.1 * xk[i];

                    if (ak[i] < Lk[i])
                    {
                        Console.WriteLine("Error: ak<Lk");
                        testbool = false;
                    }
                    if (bk[i] > Uk[i])
                    {
                        Console.WriteLine("Error: bk<Uk");
                        testbool = false;
                    }
                }
            }
            else
            {
                for (i = 0; i <= OptiSolidNumb - 1; i++)
                {
                    Lk[i] = xk[i] - 0.5;
                    Uk[i] = xk[i] + 0.5;
                    ak[i] = 0.9 * Lk[i] + 0.1 * xk[i];
                    bk[i] = 0.9 * Uk[i] + 0.1 * xk[i];

                    if (ak[i] < Lk[i])
                    {
                        Console.WriteLine("Error: ak<Lk");
                        testbool = false;
                    }
                    if (bk[i] > Uk[i])
                    {
                        Console.WriteLine("Error: bk<Uk");
                        testbool = false;
                    }
                }
            }

            //MMA calculation of pk qk
            for (i = 0; i <= OptiSolidNumb - 1; i++)
            {
                if (derfil[i] > 0)
                {
                    pk[i] = (Uk[i] - xk[i]) * (Uk[i] - xk[i]) * derfil[i];
                    qk[i] = 0;
                    Console.WriteLine("Error: positive derivative");
                    testbool = false;
                }
                else if (derfil[i] < 0)
                {
                    pk[i] = 0;
                    qk[i] = -(xk[i] - Lk[i]) * (xk[i] - Lk[i]) * derfil[i];
                }
            }
            //MMA solve langrangian duality
            al = 0.0;
            bl = Math.Pow(10, 10);
            testlamda = false;
            do
            {
                lamda = (al + bl) / 2;
                metrx = 0;
                for (i = 0; i <= OptiSolidNumb - 1; i++)
                {
                    xnew[i] = Math.Max(Math.Max(0, ak[i]), Math.Max(xk[i] - 0.2, Math.Min(Math.Min(1, bk[i]), Math.Min(xk[i] + 0.2, (Math.Sqrt(qk[i] / lamda) + Lk[i])))));
                    SHX = 0.0;
                    SH = 0.0;
                    for (j = 0; j <= 26; j++)
                    {
                        if ((Map[i, j] - 1) >= 0)
                        {
                            SHX = SHX + Dist[i, j] * xnew[Map[i, j] - 1];
                            SH = SH + Dist[i, j];
                        }
                        else if ((Map[i, j] - 1) != -1)
                        {
                            Console.WriteLine("Error: Map[i,j]-1 is {0}", Map[i, j] - 1);
                        }
                    }
                    xnewfil[i] = SHX / SH;
                    metrx = metrx + xnewfil[i];
                }
                if ((bl - al) / (al + bl) < Math.Pow(10, -3))
                {
                    testlamda = true;
                }
                else if (metrx > V)
                {
                    al = lamda;
                }
                else
                {
                    bl = lamda;
                }
            } while (testlamda == false);
        }
        public void TestProcedure(ref bool testpro, int OptiSolidNumb, int epanal, int fixedepanal, double[] xk, double[] xnew, out double max)
        {
            int i;
            max = -1.0;
            for (i = 0; i <= OptiSolidNumb - 1; i++)
            {
                if (Math.Abs(xk[i] - xnew[i]) > max)
                {
                    max = Math.Abs(xk[i] - xnew[i]);
                }
            }
            if (max < Math.Pow(10, -3) || (epanal > fixedepanal))
            {
                testpro = true;
            }
        }
        public void DeleteSolid(ref SAP2000v15.cSapModel Model, string[] OptiSolidName, int OptiSolidNumb, double[] xkfil)
        {
            int i;
            int ret;
            double temp_double;
            string temp_string;
            ret = Model.SetModelIsLocked(false);
            for (i = 0; i <= OptiSolidNumb - 1; i++)
            {
                if (xkfil[i] > 0.6)
                {
                    temp_double = Math.Round(xkfil[i] * xkfil[i] * xkfil[i], 2);
                    temp_string = "OptiSolidProp" + string.Format("{0:N2}", temp_double);
                    ret = Model.SolidObj.SetProperty(OptiSolidName[i], temp_string, 0);
                }
                else
                {
                    ret = Model.SolidObj.Delete(OptiSolidName[i], 0);
                }
            }
        }
        
        //MAIN METHOD
        static void Main(string[] args)
        {
            //DATA defined by the user
            string loadcase = "load1";    //the procedure works only for one loadcase. 
            double xmin = 0.0;
            double xmax = 2.5;
            double ymin = 0.0;
            double ymax = 0.3;
            double zmin = 0.0;
            double zmax = 1.0;
            double diakr = 0.1;
            double rmin = Math.Sqrt(3);  //the filter works for rmin <= Math.Sqrt(3)! For bigger filter check method CreateMapDist and especially the declaration of Map and Dist
            double arx = 0.3;           //volume %
            int fixedepanal = 150;       //maximum iterations

            //calculate number of elements
            int nlx = Convert.ToInt32(Math.Abs((xmax - xmin)) / diakr);
            int nly = Convert.ToInt32(Math.Abs((ymax - ymin)) / diakr);
            int nlz = Convert.ToInt32(Math.Abs((zmax - zmin)) / diakr);
            int arstoix;
            arstoix = nlx * nly * nlz;

            //INITIALIZATION
            SAP2000v15.SapObject Object;
            SAP2000v15.cSapModel Model;
            Object = new SAP2000v15.SapObject();
            Object.ApplicationStart(SAP2000v15.eUnits.kN_m_C, true, "C:\\Users\\nickvas\\Documents\\SAP\\testexample\\testexample.sdb"); //CHANGE PATH!
            Model = Object.SapModel;
            int ret;
            int i;
            bool testbool = true;

            //Create Materials and Solid properties
            Program P = new Program();
            P.CreateMatSolid(ref Model, "Europe", "EN 1992-1-1 per EN 206-1", "C30/37");

            //Read SAP Model 
            string[] PointName;
            int PointNumb;
            string[] SolidName;
            int SolidNumb;
            string[] FrameName;
            int FrameNumb;
            string[] AreaName;
            int AreaNumb;
            string[] LoadName;
            int LoadNumb;
            P.ReadModel(ref Model, out PointName, out PointNumb, out SolidName, out SolidNumb, out FrameName, out FrameNumb, out AreaName, out AreaNumb, out LoadName, out LoadNumb);

            //Create Array OptiSolidName
            int OptiSolidNumb;
            string[] OptiSolidName;
            P.CreateOptiSolidName(ref Model, ref testbool, arstoix, xmin, xmax, ymin, ymax, zmin, zmax, diakr, out OptiSolidNumb, out OptiSolidName);

            //Create Arrays Map and Dist
            int[,] Map;
            double[,] Dist;
            P.CreateMapDist(ref Model, OptiSolidNumb, nlx, nly, nlz, diakr, rmin, out Map, out Dist);

            //declaration of x
            double[] xkm2;
            xkm2 = new double[OptiSolidNumb];
            double[] xkm1;
            xkm1 = new double[OptiSolidNumb];
            double[] xk;
            xk = new double[OptiSolidNumb];
            double[] xkfil;
            xkfil = new double[OptiSolidNumb];
            double[] xnew;
            xnew = new double[OptiSolidNumb];
            double[] xnewfil;
            xnewfil = new double[OptiSolidNumb];
            //declaration of der
            double[] der;
            der = new double[OptiSolidNumb];
            double[] derfil;
            derfil = new double[OptiSolidNumb];
            double V = 0.0;
            for (i = 0; i <= OptiSolidNumb - 1; i++)
            {
                xkfil[i] = arx;
                xk[i] = arx;
                V = V + xkfil[i];
            }
            //rest of declarations
            double max;
            double comp;
            double comp2;
            
            //START OF ITERATIVE PROCCESS
            int epanal = 0;
            bool testpro = false;
            do
            {
                epanal = epanal + 1;
                Console.WriteLine("\n\n{0} iteration", epanal);

                //Update model based on array xkfil
                P.UpdateModel(ref Model, OptiSolidName, OptiSolidNumb, xkfil);

                //Set Run Case Flag
                ret = Model.Analyze.SetRunCaseFlag(loadcase, true, false);
                for (i = 0; i <= LoadNumb - 1; i++)
                {
                    if (LoadName[i] != loadcase)
                    {
                        ret = Model.Analyze.SetRunCaseFlag(LoadName[i], false, false);
                    }
                }

                //Run Analysis
                ret = Model.Analyze.RunAnalysis();

                //options for results
                ret = Model.Results.Setup.DeselectAllCasesAndCombosForOutput();
                ret = Model.Results.Setup.SetCaseSelectedForOutput(loadcase);

                if ((FrameNumb != 0) || (AreaNumb != 0) || (SolidNumb != OptiSolidNumb))
                {
                    //Calculate Compliance
                    comp = P.CalcCompliance(ref Model, ref testbool, SolidName, SolidNumb, FrameName, FrameNumb, AreaName, AreaNumb);
                }

                //Calculate Derivatives
                P.CalcDerivative(ref Model, ref testbool, OptiSolidName, OptiSolidNumb, xkfil, out der, out comp2);

                //Filter Derivatives
                derfil = P.FilterDer(ref testbool, OptiSolidNumb, Map, Dist, xkfil, der);
                
                //MMA 
                P.MMA(ref testbool, epanal, OptiSolidNumb, xk, xkm1, xkm2, derfil, Map, Dist, V, out xnew, out xnewfil);

                //Find max and testpro
                P.TestProcedure(ref testpro, OptiSolidNumb, epanal, fixedepanal, xk, xnew, out max);            
                           
                Console.WriteLine("the Compliance of all solid elements that are being optimized is: {0}", comp2);
                Console.WriteLine("no error in the procedure: {0}", testbool);
                Console.WriteLine("max is {0}", max);
                Console.WriteLine("end of an iteration");

                //MMA change values in arrays to prepare the next iteration
                for (i = 0; i <= OptiSolidNumb - 1; i++)
                {
                    xkm2[i] = xkm1[i];
                    xkm1[i] = xk[i];
                    xk[i] = xnew[i];
                    xkfil[i] = xnewfil[i];
                }

            } while (testpro == false);//END OF ITERATIVE PROCCESS

            //Delete Solids
            P.DeleteSolid(ref Model, OptiSolidName, OptiSolidNumb, xkfil);

            Console.WriteLine("The procedure is completed");
            Console.ReadKey();
        }
    }
}
