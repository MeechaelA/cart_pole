
pub struct CartPole{
    gravity: f64,

    mass_cart: f64,
    mass_pole: f64,

    length_pole: f64,

}

impl Default for CartPole{
    fn default() -> Self {
        CartPole { 
            gravity: 9.81,
            mass_cart: 0.0,
            mass_pole: 0.0, 
            length_pole: 0.0,
        }
    }
}



impl CartPole{
    pub fn new(mass_cart: f64, mass_pole: f64, length_pole: f64, gravity:f64)->Self{
        CartPole{
            mass_cart: mass_cart,
            mass_pole: mass_pole,
            length_pole: length_pole,
            gravity: gravity
        }
    }
}
//     pub fn eqn_mot_matrix_rhs(&mut self, angle: f64, angle_vel: f64,  force: f64) -> Array1<f64>{
//         let mut b:Array1<f64> = arr1(&[-self.gravity * angle.sin(), force + self.mass_pole * self.length_pole * angle_vel.powf(2.0) * angle.sin()]);
//         println!("b:{:?}", b);
//         return b;
//     }

//     pub fn eqn_mot_matrix_lhs(&mut self, angle: f64) -> Array2<f64>{
//         let eqn_mot_mat_a: Array2<f64> = arr2(&self.eqn_mot_matrix_form_lhs_a(angle));
//         return eqn_mot_mat_a;
//     }


//     pub fn eqn_mot_matrix_form_lhs_a(&mut self, angle: f64) -> [[f64; 2]; 2]{
//         let lhs_a = 
//             [
//                 [angle.cos(), self.length_pole], 
//                 [self.mass_cart+self.mass_pole , self.mass_pole * self.length_pole * angle.cos()]
//             ];
//         return lhs_a;
//     }

//     pub fn eqn_mot_matrix_form_lhs_b(&mut self, cart_accel: f64, angle_accel:f64) -> [[f64; 1]; 2]{
//         let lhs_b: [[f64; 1]; 2] = [[cart_accel], [angle_accel]];
//         return lhs_b;
//     }

//     pub fn eqn_mot_matrix_form_rhs(&mut self, force:f64, angle: f64, angle_vel: f64) -> [f64; 2]{
//         let eqn_mot_mat_rhs = 
//             [
//                 -1.0 * self.gravity * angle.cos(),
//                 force + self.mass_pole * self.length_pole * angle_vel.powf(2.0) * angle.sin(),
//             ];
//         return eqn_mot_mat_rhs
//     }
// }