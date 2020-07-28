package gauss;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import static javax.swing.JOptionPane.ERROR_MESSAGE;
import javax.swing.table.DefaultTableModel;

class MyCustomFilter extends javax.swing.filechooser.FileFilter {
    @Override
    public boolean accept(File file) {
        // Allow only directories, or files with ".txt" extension
        return file.isDirectory() || file.getAbsolutePath().endsWith(".txt");
    }
    @Override
    public String getDescription() {
        // This description will be displayed in the dialog,
        return "Text documents (*.txt)";
    }
}

// === Window class ===
public class Okno extends javax.swing.JFrame {
    
    // === Common method class (nested in class 'Okno') ===
    class SolveMethod {
        final double EPS = 0.001;                 
        
        protected boolean[] okLines;                        // for matrix degeneracy checking
        protected int dim;                                  // 'protected' - visible to the derived classes
        protected double[][] a;
        
        private String methodName;
        
        // Common part for derived classes constructors
        private void initMethod(String name) {              // 'private' - because it's used in the child' constructor (is NOT visible directly)
            dim = (int) spinDimention.getValue();           // ...in child class call it by 'super.initMethod'
            okLines = new boolean[dim];                     // autoinit with false
            methodName = name;
            
            a = new double[dim][dim+1];
            for (int i = 0; i < dim; i++)
                for (int j = 0; j < dim+1; j++) 
                    a[i][j] = (double) table.getValueAt(i,j);
        }
        
        protected void print(String s) {                    
            //System.out.print(s);                          
            textResult.append(s);                           
        }
    
        // fix -0.000
        protected String str(double x) {
            if (Math.abs(x) < EPS) 
                x = Math.abs(x);
            return String.format("%8.4f ", x);
        }
    
        protected void printMatrix(String title) {
            print(title + ":\n");
            for (int row = 0; row < dim; row++) {
                for (int col = 0; col < dim+1; col++)
                    print( str(a[row][col]) );
                print("\n");
            }
            print("\n");
        }
    
        // Partial pivoting realization
        protected boolean pivoting(int rowCur) {
            int colCur = rowCur;
            
            int rowPivot = rowCur;                          // find max element in column 
            for (int row = rowCur + 1; row < dim; row++)    // down from current
                if (Math.abs( a[row][colCur]) > Math.abs(a[rowPivot][colCur]) ) rowPivot = row;

            if (Math.abs ( a[rowPivot][colCur]) < EPS )      // if only zeros in column below (a[rowPivot][colCur] = 0)
                return false;     

            okLines[rowCur] = true;

            for (int col = 0; col < dim + 1; col++) {       // transpose rows a[rowCur] and a[rowPivot]
                double buf = a[rowCur][col];
                a[rowCur][col] = a[rowPivot][col];
                a[rowPivot][col] = buf;
            }

            return true;
        }
        
        // Check matrix for degeneracy
        protected boolean checkResults() {
            String notOneSolutionInfo = "";

            for (int row = 0; row < dim; row++) {
                if ( !okLines[row] ) {
                    if (Math.abs( a[row][dim]) > EPS ) {      // a[row][dim] !=0
                        notOneSolutionInfo = "The system is inconsistent";
                        break;
                    } else 
                        notOneSolutionInfo = "The system has a general solution";
                }
            }
            if (notOneSolutionInfo.length() == 0) return true;
                
            print(notOneSolutionInfo);
            return false;
        }
        
        void doMethod() {
            textResult.setText("");
            print("Solve the system of linear equations\n");
            print("by " + methodName + "\n\n");
            printMatrix("Source matrix");
        }
    }
    // === End of common method class ===

    // === Gauss method class ===
    class GaussMethod extends SolveMethod {
        GaussMethod() {
            super.initMethod("Gauss method");
        }
        
        @Override
        void doMethod() {
            super.doMethod();
            
            // Matrix transformation
            for (int i = 0; i < dim; i++) {
                if ( !pivoting(i) ) continue;

                double c = a[i][i];

                // Handle every row below row #i
                for (int row = i + 1; row < dim; row++) {
                    double d = a[row][i];
                    for (int col = 0; col < dim+1; col++) 
                        a[row][col] -= a[i][col] / c * d;
                }

                printMatrix( "Iteration #" + (i+1) );
            }
            // End of matrix transformation
            
            if ( !checkResults() ) return;
            
            // Calculate results
            double[] x = new double[dim+1];
            
            for(int i=dim-1; i>=0; i--){
                double s = 0;
                for(int j = i; j<dim; j++)
                    s += a[i][j] * x[j];
                x[i] = (a[i][dim] - s) / a[i][i];
            }
            
            for (int i = 0; i < dim; i++) 
                print( String.format("x%d = %s\n", i+1, str(x[i]) ) );
        }
    }
    // === End of Gauss method class ===

    // === Gauss-Jordan method class ===
    class JordanMethod extends SolveMethod {
        JordanMethod() {
            super.initMethod("Gauss-Jordan method");
        }
        
        @Override
        void doMethod() {
            super.doMethod();
            
            // Matrix transformation
            for (int i = 0; i < dim; i++) {
                if ( !pivoting(i) ) continue;

                // Divide elements of row a[rowCur] to diagonal element from this row
                double c = a[i][i];
                for (int col = 0; col < dim + 1; col++) 
                    a[i][col] = a[i][col] / c;

                // Handle every row exept row #i
                for (int row = 0; row < dim; row++)
                    if (row != i) {
                        c = a[row][i];
                        for (int col = i; col < dim+1; col++) 
                            a[row][col] -= a[i][col] * c;
                    }

                printMatrix( "Iteration #" + (i+1) );
            }
            // End of matrix transformation
        
            if ( checkResults() )
                for (int i = 0; i < dim; i++) 
                    print( String.format("x%d = %s\n", i+1, str(a[i][dim]) ) );
        }
    }    
    // === End of Gauss-Jordan method class ===
    
    // === Window class realization  ===
    Okno() {
        initComponents();
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        fileLoad = new javax.swing.JFileChooser();
        panelInput = new javax.swing.JPanel();
        labelNumber = new javax.swing.JLabel();
        spinDimention = new javax.swing.JSpinner();
        btnSolveGauss = new javax.swing.JButton();
        panelTable = new javax.swing.JPanel();
        scroolTable = new javax.swing.JScrollPane();
        table = new javax.swing.JTable();
        btnTableFromFile = new javax.swing.JButton();
        btnTableClear = new javax.swing.JButton();
        btnSolveJordan = new javax.swing.JButton();
        panelResult = new javax.swing.JPanel();
        scrollResult1 = new javax.swing.JScrollPane();
        textResult = new javax.swing.JTextArea();
        panelTitle = new javax.swing.JPanel();
        labelTitle = new javax.swing.JLabel();
        menuBar = new javax.swing.JMenuBar();
        menuFile = new javax.swing.JMenu();
        menuLoadCoefs = new javax.swing.JMenuItem();
        menuTableClear = new javax.swing.JMenuItem();
        menuSeparator1 = new javax.swing.JPopupMenu.Separator();
        menuSolveGauss = new javax.swing.JMenuItem();
        menuSolveJordan = new javax.swing.JMenuItem();
        menuSeparator2 = new javax.swing.JPopupMenu.Separator();
        menuExit = new javax.swing.JMenuItem();

        fileLoad.setDialogTitle("Select source file");
        fileLoad.setFileFilter(new MyCustomFilter());
        fileLoad.setSelectedFile(new java.io.File("D:\\Delo\\NetBeans\\03 Gauss\\Gauss\\Test11 (4).txt"));

        setDefaultCloseOperation(javax.swing.WindowConstants.DO_NOTHING_ON_CLOSE);
        setTitle("Gauss method");
        setPreferredSize(new java.awt.Dimension(800, 650));
        setResizable(false);
        addWindowListener(new java.awt.event.WindowAdapter() {
            public void windowClosing(java.awt.event.WindowEvent evt) {
                formWindowClosing(evt);
            }
            public void windowOpened(java.awt.event.WindowEvent evt) {
                formWindowOpened(evt);
            }
        });

        panelInput.setBorder(javax.swing.BorderFactory.createTitledBorder(null, " User input: ", javax.swing.border.TitledBorder.DEFAULT_JUSTIFICATION, javax.swing.border.TitledBorder.DEFAULT_POSITION, new java.awt.Font("Tahoma", 0, 11), new java.awt.Color(0, 51, 51))); // NOI18N

        labelNumber.setFont(new java.awt.Font("Tahoma", 0, 12)); // NOI18N
        labelNumber.setText("Number of equations:");

        spinDimention.setModel(new javax.swing.SpinnerNumberModel(3, 2, 9, 1));
        spinDimention.addChangeListener(new javax.swing.event.ChangeListener() {
            public void stateChanged(javax.swing.event.ChangeEvent evt) {
                spinDimentionStateChanged(evt);
            }
        });

        btnSolveGauss.setFont(new java.awt.Font("Tahoma", 1, 12)); // NOI18N
        btnSolveGauss.setText("Solve by Gauss method");
        btnSolveGauss.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnSolveGaussActionPerformed(evt);
            }
        });

        panelTable.setBorder(javax.swing.BorderFactory.createTitledBorder(null, "Source table:", javax.swing.border.TitledBorder.DEFAULT_JUSTIFICATION, javax.swing.border.TitledBorder.DEFAULT_POSITION, new java.awt.Font("Tahoma", 0, 12))); // NOI18N

        table.setModel(new javax.swing.table.DefaultTableModel(
            new Object [][] {
                {null}
            },
            new String [] {
                "Заголовок 1"
            }
        ) {
            Class[] types = new Class [] {
                java.lang.Double.class
            };

            public Class getColumnClass(int columnIndex) {
                return types [columnIndex];
            }
        });
        table.setFillsViewportHeight(true);
        scroolTable.setViewportView(table);

        btnTableFromFile.setText("Load from text file...");
        btnTableFromFile.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnTableFromFileActionPerformed(evt);
            }
        });

        btnTableClear.setText("Clear");
        btnTableClear.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnTableClearActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout panelTableLayout = new javax.swing.GroupLayout(panelTable);
        panelTable.setLayout(panelTableLayout);
        panelTableLayout.setHorizontalGroup(
            panelTableLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(scroolTable, javax.swing.GroupLayout.DEFAULT_SIZE, 305, Short.MAX_VALUE)
            .addGroup(panelTableLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(btnTableFromFile)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addComponent(btnTableClear)
                .addContainerGap())
        );
        panelTableLayout.setVerticalGroup(
            panelTableLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, panelTableLayout.createSequentialGroup()
                .addComponent(scroolTable, javax.swing.GroupLayout.PREFERRED_SIZE, 245, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addGroup(panelTableLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(btnTableFromFile)
                    .addComponent(btnTableClear)))
        );

        btnSolveJordan.setText("Solve by Gauss-Jordan method (for check)");
        btnSolveJordan.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnSolveJordanActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout panelInputLayout = new javax.swing.GroupLayout(panelInput);
        panelInput.setLayout(panelInputLayout);
        panelInputLayout.setHorizontalGroup(
            panelInputLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
            .addComponent(panelTable, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
            .addGroup(javax.swing.GroupLayout.Alignment.LEADING, panelInputLayout.createSequentialGroup()
                .addGap(16, 16, 16)
                .addComponent(labelNumber)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(spinDimention, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap(141, Short.MAX_VALUE))
            .addGroup(panelInputLayout.createSequentialGroup()
                .addGap(0, 0, Short.MAX_VALUE)
                .addGroup(panelInputLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                    .addComponent(btnSolveGauss, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(btnSolveJordan, javax.swing.GroupLayout.DEFAULT_SIZE, 267, Short.MAX_VALUE))
                .addGap(24, 24, 24))
        );
        panelInputLayout.setVerticalGroup(
            panelInputLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(panelInputLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(panelInputLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(labelNumber)
                    .addComponent(spinDimention, javax.swing.GroupLayout.PREFERRED_SIZE, 28, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(panelTable, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(25, 25, 25)
                .addComponent(btnSolveGauss, javax.swing.GroupLayout.PREFERRED_SIZE, 38, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(btnSolveJordan, javax.swing.GroupLayout.PREFERRED_SIZE, 37, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap(49, Short.MAX_VALUE))
        );

        panelResult.setBorder(javax.swing.BorderFactory.createTitledBorder(null, " Information & results: ", javax.swing.border.TitledBorder.DEFAULT_JUSTIFICATION, javax.swing.border.TitledBorder.DEFAULT_POSITION, new java.awt.Font("Tahoma", 0, 11), new java.awt.Color(0, 51, 51))); // NOI18N

        textResult.setColumns(20);
        textResult.setFont(new java.awt.Font("Consolas", 0, 12)); // NOI18N
        textResult.setRows(5);
        scrollResult1.setViewportView(textResult);

        javax.swing.GroupLayout panelResultLayout = new javax.swing.GroupLayout(panelResult);
        panelResult.setLayout(panelResultLayout);
        panelResultLayout.setHorizontalGroup(
            panelResultLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(panelResultLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(scrollResult1, javax.swing.GroupLayout.DEFAULT_SIZE, 381, Short.MAX_VALUE)
                .addContainerGap())
        );
        panelResultLayout.setVerticalGroup(
            panelResultLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(scrollResult1)
        );

        labelTitle.setFont(new java.awt.Font("Arial", 1, 18)); // NOI18N
        labelTitle.setForeground(new java.awt.Color(0, 0, 204));
        labelTitle.setText("Solve the system of linear equations by Gauss method");

        javax.swing.GroupLayout panelTitleLayout = new javax.swing.GroupLayout(panelTitle);
        panelTitle.setLayout(panelTitleLayout);
        panelTitleLayout.setHorizontalGroup(
            panelTitleLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(panelTitleLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(labelTitle)
                .addContainerGap(19, Short.MAX_VALUE))
        );
        panelTitleLayout.setVerticalGroup(
            panelTitleLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, panelTitleLayout.createSequentialGroup()
                .addGap(0, 19, Short.MAX_VALUE)
                .addComponent(labelTitle))
        );

        menuFile.setText("File");

        menuLoadCoefs.setText("Load coefs...");
        menuLoadCoefs.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                menuLoadCoefsActionPerformed(evt);
            }
        });
        menuFile.add(menuLoadCoefs);

        menuTableClear.setText("Clear matrix");
        menuTableClear.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                menuTableClearActionPerformed(evt);
            }
        });
        menuFile.add(menuTableClear);
        menuFile.add(menuSeparator1);

        menuSolveGauss.setText("Solve by Gauss method");
        menuSolveGauss.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                menuSolveGaussActionPerformed(evt);
            }
        });
        menuFile.add(menuSolveGauss);

        menuSolveJordan.setText("Solve by Gauss-Jordan method");
        menuSolveJordan.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                menuSolveJordanActionPerformed(evt);
            }
        });
        menuFile.add(menuSolveJordan);
        menuFile.add(menuSeparator2);

        menuExit.setText("Exit");
        menuExit.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                menuExitActionPerformed(evt);
            }
        });
        menuFile.add(menuExit);

        menuBar.add(menuFile);

        setJMenuBar(menuBar);

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(panelTitle, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(0, 0, Short.MAX_VALUE))
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(panelInput, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(panelResult, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addContainerGap())
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                .addComponent(panelTitle, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(18, 18, 18)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(panelResult, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(panelInput, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addContainerGap())
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents

    // Load matrix coefficients from file
    private void loadCoefsFromFile() {
        if (fileLoad.showOpenDialog(this) != JFileChooser.APPROVE_OPTION) return;
            
        int dim = (int) spinDimention.getValue();
        File file = fileLoad.getSelectedFile();
        
        try {
            redefineTable();
            
            FileReader fileReader = new FileReader( file.getAbsolutePath() );
            textResult.read( fileReader, null );
            
            String[] textLines = textResult.getText().split("\n");
            int rows = Math.min(dim,textLines.length);

            for (int i=0; i<rows; i++) {
                String textLine = textLines[i].trim().replaceAll("( |\t|;)+", " ");
                String[] coefs  = textLine.split(" ");
                int cols = Math.min(dim+1,coefs.length);
                for (int j=0; j<cols; j++) {
                    try {
                        double x = Double.parseDouble(coefs[j].trim());
                        table.setValueAt(x,i,j);
                    } catch (Exception ex){}
                }
            }
        } catch (IOException ex) {
            JOptionPane.showMessageDialog(this, "Problem accessing file:\n" + file.getAbsolutePath(), "Error opening file!", ERROR_MESSAGE);
        }
    }
    
    // Set new dimention to table & fill table by zeros
    private void redefineTable() {
        int dim = (int) spinDimention.getValue();
        
        Object[][] data = new Double[dim][dim+1];
        Object[] headers = new String [dim+1];
        Class[] typesDouble = new Class [dim+1];
        
        for (int j = 0; j < dim+1; j++) {
            headers[j] = "";
            typesDouble[j] = java.lang.Double.class;
        }
        
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim+1; j++)
                data[i][j] = 0.0;
        
        DefaultTableModel newModel = new DefaultTableModel(data, headers) {
            Class[] types = typesDouble;

            @Override
            public Class getColumnClass(int columnIndex) {
                return types [columnIndex];
            }
        };
        
        table.setModel(newModel);
        textResult.setText("");
    }
    
    // Exit from program
    private void closeProgram(java.awt.AWTEvent evt) {
        String message = "Are you really want to exit?";
        String title = "Exit";
        int reply = JOptionPane.showConfirmDialog(this, message, title, JOptionPane.YES_NO_OPTION);
        if (reply == JOptionPane.YES_OPTION) {
            System.exit(0);
        }
    }
    
    private void btnSolveGaussActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnSolveGaussActionPerformed
        new GaussMethod().doMethod();
    }//GEN-LAST:event_btnSolveGaussActionPerformed
    
    private void btnTableFromFileActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnTableFromFileActionPerformed
        this.loadCoefsFromFile();
    }//GEN-LAST:event_btnTableFromFileActionPerformed

    private void menuLoadCoefsActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_menuLoadCoefsActionPerformed
        this.loadCoefsFromFile();
    }//GEN-LAST:event_menuLoadCoefsActionPerformed

    private void menuSolveGaussActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_menuSolveGaussActionPerformed
        new GaussMethod().doMethod();
    }//GEN-LAST:event_menuSolveGaussActionPerformed

    private void menuExitActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_menuExitActionPerformed
        this.closeProgram(evt);
    }//GEN-LAST:event_menuExitActionPerformed

    private void formWindowClosing(java.awt.event.WindowEvent evt) {//GEN-FIRST:event_formWindowClosing
        this.closeProgram(evt);
    }//GEN-LAST:event_formWindowClosing

    private void formWindowOpened(java.awt.event.WindowEvent evt) {//GEN-FIRST:event_formWindowOpened
        this.redefineTable();
    }//GEN-LAST:event_formWindowOpened

    private void spinDimentionStateChanged(javax.swing.event.ChangeEvent evt) {//GEN-FIRST:event_spinDimentionStateChanged
        this.redefineTable();
    }//GEN-LAST:event_spinDimentionStateChanged

    private void btnTableClearActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnTableClearActionPerformed
        this.redefineTable();
    }//GEN-LAST:event_btnTableClearActionPerformed

    private void btnSolveJordanActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnSolveJordanActionPerformed
        new JordanMethod().doMethod();
    }//GEN-LAST:event_btnSolveJordanActionPerformed

    private void menuSolveJordanActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_menuSolveJordanActionPerformed
        new JordanMethod().doMethod();
    }//GEN-LAST:event_menuSolveJordanActionPerformed

    private void menuTableClearActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_menuTableClearActionPerformed
        this.redefineTable();
    }//GEN-LAST:event_menuTableClearActionPerformed
    
    // <editor-fold defaultstate="collapsed" desc="Variables declaration">
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton btnSolveGauss;
    private javax.swing.JButton btnSolveJordan;
    private javax.swing.JButton btnTableClear;
    private javax.swing.JButton btnTableFromFile;
    private javax.swing.JFileChooser fileLoad;
    private javax.swing.JLabel labelNumber;
    private javax.swing.JLabel labelTitle;
    private javax.swing.JMenuBar menuBar;
    private javax.swing.JMenuItem menuExit;
    private javax.swing.JMenu menuFile;
    private javax.swing.JMenuItem menuLoadCoefs;
    private javax.swing.JPopupMenu.Separator menuSeparator1;
    private javax.swing.JPopupMenu.Separator menuSeparator2;
    private javax.swing.JMenuItem menuSolveGauss;
    private javax.swing.JMenuItem menuSolveJordan;
    private javax.swing.JMenuItem menuTableClear;
    private javax.swing.JPanel panelInput;
    private javax.swing.JPanel panelResult;
    private javax.swing.JPanel panelTable;
    private javax.swing.JPanel panelTitle;
    private javax.swing.JScrollPane scrollResult1;
    private javax.swing.JScrollPane scroolTable;
    private javax.swing.JSpinner spinDimention;
    private javax.swing.JTable table;
    private javax.swing.JTextArea textResult;
    // End of variables declaration//GEN-END:variables
    // </editor-fold>                        
}

