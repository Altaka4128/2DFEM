using MathNet.Numerics.LinearAlgebra.Double;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace _2DFEM
{   
    class Shell1dElement : FEM_Element
    {
        const int NodeNum = 4;
        const int IntegralPoints = 4;

        private double Young;       // ヤング率
        private double Poisson;     // ポアソン比
        private double Thickness;   // 厚み

        private double[] Xi = new double[IntegralPoints];    // ξ(定数)
        private double[] Et = new double[IntegralPoints];    // η(定数)
        private double[] w_i = new double[IntegralPoints];   // 積分点重みwi
        private double[] w_j = new double[IntegralPoints];   // 積分点重みwj

        private DenseMatrix DMatrix;                                       // Dマトリックス
        private DenseMatrix[] BMatrix = new DenseMatrix[IntegralPoints];   // Bマトリックス
        private DenseMatrix[] JMatrix = new DenseMatrix[IntegralPoints];   // ヤコビ行列
        private DenseMatrix KeMatrix;                                      // Keマトリックス 
        public DenseVector[] IntegralStrainVector   // ひずみベクトル
        {
            get;
            private set;
        }
        public DenseVector[] IntegralStressVector   // 応力ベクトル
        {
            get;
            private set;
        }
        public DenseVector AveStressVector   // 要素内の平均応力ベクトル
        {
            get;
            private set;
        }
        public override DenseVector[] NodeStressVector
        {
            get;
            protected set;
        }

        public Shell1dElement(Node[] nodes,
                              double thickness,
                              double young,
                              double poisson)
        {

            if (nodes.Length != NodeNum)
            {
                return;
            }

            Nodes = nodes;
            Thickness = thickness;
            Young = young;
            Poisson = poisson;
            BMatrix = new DenseMatrix[IntegralPoints];
            IntegralStrainVector = new DenseVector[IntegralPoints];
            IntegralStressVector = new DenseVector[IntegralPoints];
            NodeStressVector = new DenseVector[NodeNum];

            // ξ(定数)を初期化する
            Xi[0] = -1 / Math.Sqrt(3);
            Xi[1] = 1 / Math.Sqrt(3);
            Xi[2] = -1 / Math.Sqrt(3);
            Xi[3] = 1 / Math.Sqrt(3);

            // η(定数)を初期化する
            Et[0] = -1 / Math.Sqrt(3);
            Et[1] = -1 / Math.Sqrt(3);
            Et[2] = 1 / Math.Sqrt(3);
            Et[3] = 1 / Math.Sqrt(3);

            // 積分点重みを初期化する
            w_i[0] = 1.0;
            w_i[1] = 1.0;
            w_i[2] = 1.0;
            w_i[3] = 1.0;
            w_j[0] = 1.0;
            w_j[1] = 1.0;
            w_j[2] = 1.0;
            w_j[3] = 1.0;
        }

        // Dマトリックスを計算する
        private void makeDMatirx()
        {
            // 例外処理
            if (Young == 0.0)
            {
                return;
            }

            double coef = Young / ((1 - 2 * Poisson) * (1 + Poisson));
            double[,] dmatrixArray = new double[3, 3];
            dmatrixArray[0, 0] = 1 - Poisson;
            dmatrixArray[0, 1] = Poisson;
            dmatrixArray[0, 2] = 0;
            dmatrixArray[1, 0] = Poisson;
            dmatrixArray[1, 1] = 1 - Poisson;
            dmatrixArray[1, 2] = 0;
            dmatrixArray[2, 0] = 0;
            dmatrixArray[2, 1] = 0;
            dmatrixArray[2, 2] = (1 - 2 * Poisson) / 2;
            DMatrix = coef * DenseMatrix.OfArray(dmatrixArray);

            Console.WriteLine("Dマトリックス");
            Console.WriteLine(DMatrix);
        }

        // Bマトリックスを計算する
        private void makeBMatirx()
        {
            // 例外処理
            if (Thickness <= 0)
            {
                return;
            }

            for(int i = 0; i < IntegralPoints; i++)
            {
                // dNi/dξを計算する
                double dN1dXi = -(1 - Et[i]) / 4;
                double dN2dXi =  (1 - Et[i]) / 4;
                double dN3dXi =  (1 + Et[i]) / 4;
                double dN4dXi = -(1 + Et[i]) / 4;
                // dNi/dηを計算する
                double dN1dEt = -(1 - Xi[i]) / 4;
                double dN2dEt = -(1 + Xi[i]) / 4;
                double dN3dEt =  (1 + Xi[i]) / 4;
                double dN4dEt =  (1 - Xi[i]) / 4;
                // dx/dξ、dy/dξを計算する
                double dxdXi = dN1dXi * Nodes[0].Point.X + dN2dXi * Nodes[1].Point.X +
                               dN3dXi * Nodes[2].Point.X + dN4dXi * Nodes[3].Point.X;
                double dydXi = dN1dXi * Nodes[0].Point.Y + dN2dXi * Nodes[1].Point.Y +
                               dN3dXi * Nodes[2].Point.Y + dN4dXi * Nodes[3].Point.Y;
                // dx/dη、dy/dηを計算する
                double dxdEt = dN1dEt * Nodes[0].Point.X + dN2dEt * Nodes[1].Point.X +
                               dN3dEt * Nodes[2].Point.X + dN4dEt * Nodes[3].Point.X;
                double dydEt = dN1dEt * Nodes[0].Point.Y + dN2dEt * Nodes[1].Point.Y +
                               dN3dEt * Nodes[2].Point.Y + dN4dEt * Nodes[3].Point.Y;
                // ヤコビ行列Jを計算する
                double[,] jmatrixArray = new double[2, 2];
                jmatrixArray[0, 0] = dxdXi;
                jmatrixArray[0, 1] = dydXi;
                jmatrixArray[1, 0] = dxdEt;
                jmatrixArray[1, 1] = dydEt;
                JMatrix[i] = DenseMatrix.OfArray(jmatrixArray);
                DenseMatrix JMatrix_inv = (DenseMatrix)JMatrix[i].Inverse();
                // dNi/dx、dNi/dyを計算する
                DenseVector dN1dXi_EtVector = DenseVector.OfArray(new double[] { dN1dXi, dN1dEt });
                DenseVector dN2dXi_EtVector = DenseVector.OfArray(new double[] { dN2dXi, dN2dEt });
                DenseVector dN3dXi_EtVector = DenseVector.OfArray(new double[] { dN3dXi, dN3dEt });
                DenseVector dN4dXi_EtVector = DenseVector.OfArray(new double[] { dN4dXi, dN4dEt });
                DenseVector dN1dx_yVector = (DenseVector)JMatrix_inv.Multiply(dN1dXi_EtVector);
                DenseVector dN2dx_yVector = (DenseVector)JMatrix_inv.Multiply(dN2dXi_EtVector);
                DenseVector dN3dx_yVector = (DenseVector)JMatrix_inv.Multiply(dN3dXi_EtVector);
                DenseVector dN4dx_yVector = (DenseVector)JMatrix_inv.Multiply(dN4dXi_EtVector);

                // Bマトリックスを計算する
                double[,] bmatrixArray = new double[3, NodeNum * 2];
                bmatrixArray[0, 0] = dN1dx_yVector[0];
                bmatrixArray[0, 1] = 0;
                bmatrixArray[0, 2] = dN2dx_yVector[0];
                bmatrixArray[0, 3] = 0;
                bmatrixArray[0, 4] = dN3dx_yVector[0];
                bmatrixArray[0, 5] = 0;
                bmatrixArray[0, 6] = dN4dx_yVector[0];
                bmatrixArray[0, 7] = 0;
                bmatrixArray[1, 0] = 0;
                bmatrixArray[1, 1] = dN1dx_yVector[1];
                bmatrixArray[1, 2] = 0;
                bmatrixArray[1, 3] = dN2dx_yVector[1];
                bmatrixArray[1, 4] = 0;
                bmatrixArray[1, 5] = dN3dx_yVector[1];
                bmatrixArray[1, 6] = 0;
                bmatrixArray[1, 7] = dN4dx_yVector[1];
                bmatrixArray[2, 0] = dN1dx_yVector[1];
                bmatrixArray[2, 1] = dN1dx_yVector[0];
                bmatrixArray[2, 2] = dN2dx_yVector[1];
                bmatrixArray[2, 3] = dN2dx_yVector[0];
                bmatrixArray[2, 4] = dN3dx_yVector[1];
                bmatrixArray[2, 5] = dN3dx_yVector[0];
                bmatrixArray[2, 6] = dN4dx_yVector[1];
                bmatrixArray[2, 7] = dN4dx_yVector[0];
                BMatrix[i] = DenseMatrix.OfArray(bmatrixArray);

                Console.WriteLine("積分点" + (i + 1) + " Bマトリックス");
                Console.WriteLine(BMatrix[i]);
            }
        }

        // Keマトリックスを計算する
        public override DenseMatrix makeKeMatrix()
        {
            if(Thickness <= 0)
            {
                return null;
            }

            // Dマトリックスを計算する
            makeDMatirx();

            // Bマトリックスを計算する(Jマトリックスも途中過程で計算する)
            makeBMatirx();

            DenseMatrix KeMatrix = DenseMatrix.Create(NodeNum * 2, NodeNum * 2, 0.0);
            for (int i = 0; i < 4; i++)
            {

                double coef = Thickness * w_i[i] * w_j[i] * JMatrix[i].Determinant();
                DenseMatrix BtMatrix = (DenseMatrix)BMatrix[i].Transpose();
                DenseMatrix DtMatrix = (DenseMatrix)DMatrix.Transpose();
                DenseMatrix KepMatrix = coef * (DenseMatrix)(BtMatrix.Multiply(DtMatrix)).Multiply(BMatrix[i]);
                KeMatrix += KepMatrix;
            }

            Console.WriteLine("Keマトリックス");
            Console.WriteLine(KeMatrix);

            return KeMatrix;
        }

        // ひずみベクトルを計算する
        public override void makeStrainVector(DenseVector dispvector)
        {
            for (int i = 0; i < IntegralPoints; i++)
            {
                IntegralStrainVector[i] = (DenseVector)BMatrix[i].Multiply(dispvector);
                Console.WriteLine("積分点" + (i + 1) + "ひずみベクトル");
                Console.WriteLine(IntegralStrainVector[i]);
            }
        }

        // 応力ベクトルを計算する
        public override void makeStressVector()
        {
            for (int i = 0; i < IntegralPoints; i++)
            {
                // 例外処理
                if (IntegralStrainVector[i] == null)
                {
                    return;
                }

                IntegralStressVector[i] = (DenseVector)DMatrix.Multiply(IntegralStrainVector[i]);
                Console.WriteLine("積分点" + (i + 1) + "応力ベクトル");
                Console.WriteLine(IntegralStressVector[i]);
            }
        }


        public void makeAveStress()
        {
            // 各積分点の応力の合計値を計算する
            DenseVector sumStressVector = DenseVector.Create(3, 0.0);
            for (int i = 0; i < IntegralPoints; i++)
            {
                // 例外処理
                if (IntegralStressVector[i] == null)
                {
                    return;
                }
                sumStressVector += IntegralStressVector[i];
            }

            // 平均値を計算する
            AveStressVector = sumStressVector / (double)IntegralPoints;

            Console.WriteLine("平均応力ベクトル");
            Console.WriteLine(AveStressVector);
        }

        // 線形外挿法で節点応力を計算する
        public override void makeNodeStressVector()
        {
            // 例外処理
            for (int i = 0; i < NodeNum; i++)
            {
                if (IntegralStressVector[i] == null)
                {
                    return;
                }
            }

            // 要素の平均応力を計算する
            makeAveStress();

            // 各節点の応力を計算する(計算方法が合っているか要検討)
            double cof = Math.Sqrt(3);
            NodeStressVector[0] = AveStressVector + cof * (IntegralStressVector[0] - AveStressVector);
            NodeStressVector[1] = AveStressVector + cof * (IntegralStressVector[1] - AveStressVector);
            NodeStressVector[2] = AveStressVector + cof * (IntegralStressVector[3] - AveStressVector);
            NodeStressVector[3] = AveStressVector + cof * (IntegralStressVector[2] - AveStressVector);

            for(int i = 0; i < NodeNum; i++)
            {
                Console.WriteLine("ノード番号" + (i + 1) + "節点応力");
                Console.WriteLine(NodeStressVector[i]);
            }
        }
    }
}
