using MathNet.Numerics.LinearAlgebra.Double;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace _2DFEM
{
    class Shell2dElement : FEM_Element
    {
        const int NodeNum = 8;
        const int IntegralPoints = 9;

        private double Thickness;
        private double Young;
        private double Poisson;

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
        public override DenseVector[] NodeStressVector
        {
            get;
            protected set;
        }


        public Shell2dElement(Node[] nodes, double thickness, double young, double poisson)
        {

            if (nodes.Length != NodeNum)
            {
                return;
            }

            Nodes = nodes;
            Thickness = thickness;
            Young = young;
            Poisson = poisson;
            IntegralStrainVector = new DenseVector[IntegralPoints];
            IntegralStressVector = new DenseVector[IntegralPoints];
            NodeStressVector = new DenseVector[NodeNum];

            // ξ(定数)、η(定数)を初期化する
            double cof = Math.Sqrt(3.0 / 5.0);
            Xi[0] = -cof;
            Xi[1] = 0;
            Xi[2] = cof;
            Xi[3] = -cof;
            Xi[4] = 0;
            Xi[5] = cof;
            Xi[6] = -cof;
            Xi[7] = 0;
            Xi[8] = cof;

            Et[0] = -cof;
            Et[1] = -cof;
            Et[2] = -cof;
            Et[3] = 0;
            Et[4] = 0;
            Et[5] = 0;
            Et[6] = cof;
            Et[7] = cof;
            Et[8] = cof;

            // 積分点重みを初期化する
            double we = 5.0 / 9.0;
            double wc = 8.0 / 9.0;
            w_i[0] = we;
            w_i[1] = wc;
            w_i[2] = we;
            w_i[3] = we;
            w_i[4] = wc;
            w_i[5] = we;
            w_i[6] = we;
            w_i[7] = wc;
            w_i[8] = we;

            w_j[0] = we;
            w_j[1] = we;
            w_j[2] = we;
            w_j[3] = wc;
            w_j[4] = wc;
            w_j[5] = wc;
            w_j[6] = we;
            w_j[7] = we;
            w_j[8] = we;

            makeBMatirx();
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

            for (int i = 0; i < IntegralPoints; i++)
            {
                // dNi/dξを計算する
                double dN1dXi = (-(1 - Et[i]) * (-1 - Xi[i] - Et[i]) - (1 - Xi[i]) * (1 - Et[i])) / 4.0;
                double dN2dXi = ((1 - Et[i]) * (-1 + Xi[i] - Et[i]) + (1 + Xi[i]) * (1 - Et[i])) / 4.0;
                double dN3dXi = ((1 + Et[i]) * (-1 + Xi[i] + Et[i]) + (1 + Xi[i]) * (1 + Et[i])) / 4.0;
                double dN4dXi = (-(1 + Et[i]) * (-1 - Xi[i] + Et[i]) - (1 - Xi[i]) * (1 + Et[i])) / 4.0;
                double dN5dXi = (-2 * Xi[i] * (1 - Et[i])) / 2.0;
                double dN6dXi = (1 - Et[i] * Et[i]) / 2.0;
                double dN7dXi = (-2 * Xi[i] * (1 + Et[i])) / 2.0;
                double dN8dXi = (-(1 - Et[i] * Et[i])) / 2.0;
                // dNi/dηを計算する
                double dN1dEt = (-(1 - Xi[i]) * (-1 - Xi[i] - Et[i]) - (1 - Xi[i]) * (1 - Et[i])) / 4.0;
                double dN2dEt = (-(1 + Xi[i]) * (-1 + Xi[i] - Et[i]) - (1 + Xi[i]) * (1 - Et[i])) / 4.0;
                double dN3dEt = ((1 + Xi[i]) * (-1 + Xi[i] + Et[i]) + (1 + Xi[i]) * (1 + Et[i])) / 4.0;
                double dN4dEt = ((1 - Xi[i]) * (-1 - Xi[i] + Et[i]) + (1 - Xi[i]) * (1 + Et[i])) / 4.0;
                double dN5dEt = (-(1 - Xi[i] * Xi[i])) / 2.0;
                double dN6dEt = (-2 * Et[i] * (1 + Xi[i])) / 2.0;
                double dN7dEt = ((1 - Xi[i] * Xi[i])) / 2.0;
                double dN8dEt = (-2 * Et[i] * (1 - Xi[i])) / 2.0;
                // dx/dξ、dy/dξを計算する
                double dxdXi = dN1dXi * Nodes[0].Point.X + dN2dXi * Nodes[1].Point.X +
                               dN3dXi * Nodes[2].Point.X + dN4dXi * Nodes[3].Point.X +
                               dN5dXi * Nodes[4].Point.X + dN6dXi * Nodes[5].Point.X +
                               dN7dXi * Nodes[6].Point.X + dN8dXi * Nodes[7].Point.X;

                double dydXi = dN1dXi * Nodes[0].Point.Y + dN2dXi * Nodes[1].Point.Y +
                               dN3dXi * Nodes[2].Point.Y + dN4dXi * Nodes[3].Point.Y +
                               dN5dXi * Nodes[4].Point.Y + dN6dXi * Nodes[5].Point.Y +
                               dN7dXi * Nodes[6].Point.Y + dN8dXi * Nodes[7].Point.Y;
                // dx/dη、dy/dηを計算する
                double dxdEt = dN1dEt * Nodes[0].Point.X + dN2dEt * Nodes[1].Point.X +
                               dN3dEt * Nodes[2].Point.X + dN4dEt * Nodes[3].Point.X +
                               dN5dEt * Nodes[4].Point.X + dN6dEt * Nodes[5].Point.X +
                               dN7dEt * Nodes[6].Point.X + dN8dEt * Nodes[7].Point.X;
                double dydEt = dN1dEt * Nodes[0].Point.Y + dN2dEt * Nodes[1].Point.Y +
                               dN3dEt * Nodes[2].Point.Y + dN4dEt * Nodes[3].Point.Y +
                               dN5dEt * Nodes[4].Point.Y + dN6dEt * Nodes[5].Point.Y +
                               dN7dEt * Nodes[6].Point.Y + dN8dEt * Nodes[7].Point.Y;

                // ヤコビ行列Jを計算する
                double[,] jmatrixArray = new double[2, 2];
                jmatrixArray[0, 0] = dxdXi;
                jmatrixArray[0, 1] = dydXi;
                jmatrixArray[1, 0] = dxdEt;
                jmatrixArray[1, 1] = dydEt;
                JMatrix[i] = DenseMatrix.OfArray(jmatrixArray);
                DenseMatrix JMatrix_inv = (DenseMatrix)JMatrix[i].Inverse();
                //dNi / dx、dNi / dyを計算する
                DenseVector dN1dXi_EtVector = DenseVector.OfArray(new double[] { dN1dXi, dN1dEt });
                DenseVector dN2dXi_EtVector = DenseVector.OfArray(new double[] { dN2dXi, dN2dEt });
                DenseVector dN3dXi_EtVector = DenseVector.OfArray(new double[] { dN3dXi, dN3dEt });
                DenseVector dN4dXi_EtVector = DenseVector.OfArray(new double[] { dN4dXi, dN4dEt });
                DenseVector dN5dXi_EtVector = DenseVector.OfArray(new double[] { dN5dXi, dN5dEt });
                DenseVector dN6dXi_EtVector = DenseVector.OfArray(new double[] { dN6dXi, dN6dEt });
                DenseVector dN7dXi_EtVector = DenseVector.OfArray(new double[] { dN7dXi, dN7dEt });
                DenseVector dN8dXi_EtVector = DenseVector.OfArray(new double[] { dN8dXi, dN8dEt });
                DenseVector dN1dx_yVector = (DenseVector)JMatrix_inv.Multiply(dN1dXi_EtVector);
                DenseVector dN2dx_yVector = (DenseVector)JMatrix_inv.Multiply(dN2dXi_EtVector);
                DenseVector dN3dx_yVector = (DenseVector)JMatrix_inv.Multiply(dN3dXi_EtVector);
                DenseVector dN4dx_yVector = (DenseVector)JMatrix_inv.Multiply(dN4dXi_EtVector);
                DenseVector dN5dx_yVector = (DenseVector)JMatrix_inv.Multiply(dN5dXi_EtVector);
                DenseVector dN6dx_yVector = (DenseVector)JMatrix_inv.Multiply(dN6dXi_EtVector);
                DenseVector dN7dx_yVector = (DenseVector)JMatrix_inv.Multiply(dN7dXi_EtVector);
                DenseVector dN8dx_yVector = (DenseVector)JMatrix_inv.Multiply(dN8dXi_EtVector);

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
                bmatrixArray[0, 8] = dN5dx_yVector[0];
                bmatrixArray[0, 9] = 0;
                bmatrixArray[0, 10] = dN6dx_yVector[0];
                bmatrixArray[0, 11] = 0;
                bmatrixArray[0, 12] = dN7dx_yVector[0];
                bmatrixArray[0, 13] = 0;
                bmatrixArray[0, 14] = dN8dx_yVector[0];
                bmatrixArray[0, 15] = 0;
                bmatrixArray[1, 0] = 0;
                bmatrixArray[1, 1] = dN1dx_yVector[1];
                bmatrixArray[1, 2] = 0;
                bmatrixArray[1, 3] = dN2dx_yVector[1];
                bmatrixArray[1, 4] = 0;
                bmatrixArray[1, 5] = dN3dx_yVector[1];
                bmatrixArray[1, 6] = 0;
                bmatrixArray[1, 7] = dN4dx_yVector[1];
                bmatrixArray[1, 8] = 0;
                bmatrixArray[1, 9] = dN5dx_yVector[1];
                bmatrixArray[1, 10] = 0;
                bmatrixArray[1, 11] = dN6dx_yVector[1];
                bmatrixArray[1, 12] = 0;
                bmatrixArray[1, 13] = dN7dx_yVector[1];
                bmatrixArray[1, 14] = 0;
                bmatrixArray[1, 15] = dN8dx_yVector[1];
                bmatrixArray[2, 0] = dN1dx_yVector[1];
                bmatrixArray[2, 1] = dN1dx_yVector[0];
                bmatrixArray[2, 2] = dN2dx_yVector[1];
                bmatrixArray[2, 3] = dN2dx_yVector[0];
                bmatrixArray[2, 4] = dN3dx_yVector[1];
                bmatrixArray[2, 5] = dN3dx_yVector[0];
                bmatrixArray[2, 6] = dN4dx_yVector[1];
                bmatrixArray[2, 7] = dN4dx_yVector[0];
                bmatrixArray[2, 8] = dN5dx_yVector[1];
                bmatrixArray[2, 9] = dN5dx_yVector[0];
                bmatrixArray[2, 10] = dN6dx_yVector[1];
                bmatrixArray[2, 11] = dN6dx_yVector[0];
                bmatrixArray[2, 12] = dN7dx_yVector[1];
                bmatrixArray[2, 13] = dN7dx_yVector[0];
                bmatrixArray[2, 14] = dN8dx_yVector[1];
                bmatrixArray[2, 15] = dN8dx_yVector[0];
                BMatrix[i] = DenseMatrix.OfArray(bmatrixArray);

                Console.WriteLine("積分点" + (i + 1) + " Bマトリックス");
                Console.WriteLine(BMatrix[i]);
            }
        }

        // Keマトリックスを計算する
        public override DenseMatrix makeKeMatrix()
        {
            if (Thickness <= 0)
            {
                return null;
            }

            // Dマトリックスを計算する
            makeDMatirx();

            // Bマトリックスを計算する(Jマトリックスも途中過程で計算する)
            makeBMatirx();

            DenseMatrix KeMatrix = DenseMatrix.Create(NodeNum * 2, NodeNum * 2, 0.0);
            for (int i = 0; i < IntegralPoints; i++)
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

        public override void makeStrainVector(DenseVector dispvector)
        {
            for (int i = 0; i < IntegralPoints; i++)
            {
                IntegralStrainVector[i] = (DenseVector)BMatrix[i].Multiply(dispvector);
                Console.WriteLine("積分点" + (i + 1) + "ひずみベクトル");
                Console.WriteLine(IntegralStrainVector[i]);
            }
        }

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

        public DenseVector makeAveStress()
        {
            // 各積分点の応力の合計値を計算する
            DenseVector sumStressVector = DenseVector.Create(3, 0.0);
            for (int i = 0; i < IntegralPoints; i++)
            {
                // 例外処理
                if (IntegralStressVector[i] == null)
                {
                    return null;
                }
                sumStressVector += IntegralStressVector[i];
            }

            // 平均値を計算する
            DenseVector aveStressVector = sumStressVector / (double)IntegralPoints;

            Console.WriteLine("平均応力ベクトル");
            Console.WriteLine(aveStressVector);

            return aveStressVector;
        }

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
            DenseVector aveStressVector = makeAveStress();

            // 各節点の応力を計算する(計算方法が合っているか要検討)
            double cof = Math.Sqrt(5.0 / 3.0);
            NodeStressVector[0] = aveStressVector + cof * (IntegralStressVector[0] - aveStressVector);
            NodeStressVector[1] = aveStressVector + cof * (IntegralStressVector[2] - aveStressVector);
            NodeStressVector[2] = aveStressVector + cof * (IntegralStressVector[8] - aveStressVector);
            NodeStressVector[3] = aveStressVector + cof * (IntegralStressVector[6] - aveStressVector);
            NodeStressVector[4] = aveStressVector + cof * (IntegralStressVector[1] - aveStressVector);
            NodeStressVector[5] = aveStressVector + cof * (IntegralStressVector[5] - aveStressVector);
            NodeStressVector[6] = aveStressVector + cof * (IntegralStressVector[7] - aveStressVector);
            NodeStressVector[7] = aveStressVector + cof * (IntegralStressVector[3] - aveStressVector);

            for (int i = 0; i < NodeNum; i++)
            {
                Console.WriteLine("ノード番号" + (i + 1) + "節点応力");
                Console.WriteLine(NodeStressVector[i]);
            }
        }
    }
}
