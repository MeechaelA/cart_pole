pub enum IntegrationMessage{
    ValidStep
}

pub struct TrapezoidalIntegrator{

}

impl Default for TrapezoidalIntegrator{
    fn default() -> Self {
        TrapezoidalIntegrator {  
        }
    }
}

impl TrapezoidalIntegrator{
    pub fn integrate(&mut self, x: Vec<f64>, y_func: &dyn Fn(f64) -> f64)->f64{
        
        let mut x_0: f64 = 0.0; // x-1
        let mut x_1: f64 = 0.0; // x
        let mut d_x: f64 = 0.0;
        let mut sum: f64 = 0.0;
        for (i, value) in x.iter().enumerate(){
            if i != 0{
                x_1 = *value;
                d_x = x_1 - x_0;
                let f_x_0 = y_func(x_0);
                let f_x_1 = y_func(x_1);
                sum += (f_x_1 + f_x_0)/2.0 * d_x;
                x_0 = *value;
            }
            else{
                x_0 = *value;
            }
        }
        return sum;
    }
}

pub struct TrapezoidalIntegratorStep{
    x: Vec<f64>,
    y: Vec<f64>
}

impl Default for TrapezoidalIntegratorStep{
    fn default() -> Self {
        TrapezoidalIntegratorStep {  
            x: vec![],
            y: vec![]
        }
    }
}

impl TrapezoidalIntegratorStep{
    pub fn initialize(&mut self, x_initial: f64, y_initial: f64){
        self.x.push(x_initial);
        self.y.push(y_initial)
    }

    pub fn integrate(&mut self, new_x: f64, new_y: f64)->Result<IntegrationMessage, &'static str>{
        let mut x_0: f64 = 0.0; // x-1
        let mut d_x: f64 = 0.0;
        let mut sum: f64 = 0.0;

        let x_last = self.x.last();
        let y_last = self.y.last();
        match x_last{
            Some(x_last) =>{
                match y_last{
                    Some(y_last)=>{
                        let x_last_save = x_last.clone();
                        let y_last_save = y_last.clone();
                        self.x.push(x_last_save);
        
                        d_x = new_x - x_last_save;
                        sum += (new_y + y_last)/2.0 * d_x;
                        self.y.push(sum);
                        return Ok(IntegrationMessage::ValidStep);                      
                    }
                    None=>{
                        return Err("Check if Y Vector is greater than size 0");
                    }
                }

            }
            None=>{
                return Err("Check if X Vector is greater than size 0");
            }
        }
    }

    pub fn get_y_at(&mut self, index: usize) -> f64{
        return self.y[index];
    }
}