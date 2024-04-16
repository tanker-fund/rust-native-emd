#[macro_use]
extern crate lazy_static;
use std::sync::Mutex;
use std::rc::Rc;
use std::cell::RefCell;

const EPSILON: f64 = 1e-6;

lazy_static! {
    static ref CURRENT_NODE_ID: Mutex<usize> = Mutex::new(0);
}


#[derive(Debug, PartialEq)]
enum Color {
    WHITE,
    GRAY,
    BLACK,
}

#[derive(Debug)]
struct Basic {
    id: usize,
    row: usize,
    col: usize,
    flow: f64,
    adjacency: Vec<Rc<RefCell<Basic>>>,
    cursor: Option<usize>,
    back: Option<Rc<RefCell<Basic>>>,
    color: Color,
}
type BasicRef = Rc<RefCell<Basic>>;
type Basis = Vec<BasicRef>;

fn initialize_flow(weight_x: &Vec<f64>, weight_y: &Vec<f64>, cost: &Vec<Vec<f64>>) -> Basis {
    let n_x = weight_x.len();
    let n_y = weight_y.len();
    let mut basis: Basis = Vec::new();
    let mut remaining_x = weight_x.clone();
    let mut remaining_y = weight_y.clone();
    let mut fx: usize = 0;
    let mut fy: usize = 0;
    let mut b: usize = 0;
    let b_max: usize = n_x + n_y - 1;

    while b < b_max {
        if fx == (n_x - 1) {
            for fy in fy..n_y {
                let basic = init_basic(fx, fy, remaining_y[fy]);
                insert_basic(&mut basis, b, &basic);
                b += 1;
            }
            break;
        }
        if fy == (n_y - 1) {
            for fx in fx..n_x {
                let basic = init_basic(fx, fy, remaining_x[fx]);
                insert_basic(&mut basis, b, &basic);
                b += 1;
            }
            break;
        }
        if remaining_x[fx] <= remaining_y[fy] {
            let basic = init_basic(fx, fy, remaining_x[fx]);
            insert_basic(&mut basis, b, &basic);
            b += 1;
            remaining_y[fy] -= remaining_x[fx];
            fx += 1;
        } else {
            let basic = init_basic(fx, fy, remaining_y[fy]);
            insert_basic(&mut basis, b, &basic);
            b += 1;
            remaining_x[fx] -= remaining_y[fy];
            fy += 1;
        }
    }

    basis
}

fn init_basic(row: usize, col: usize, flow: f64) -> BasicRef {
    let mut id_guard = CURRENT_NODE_ID.lock().unwrap();
    *id_guard += 1;
    Rc::new(RefCell::new(Basic {
        id: *id_guard,
        row,
        col,
        flow,
        adjacency: Vec::<Rc<RefCell<Basic>>>::new(),
        cursor: None,
        back: None,
        color: Color::WHITE,
    }))
}

fn insert_basic(basis: &mut Basis, index: usize, node_ref: &BasicRef) {
    if basis.len() < index {
        panic!("Index is to large");
    } else if basis.len() <= index {
        basis.push(Rc::clone(node_ref));
    } else {
        basis[index] = Rc::clone(node_ref);
    }
    let mut node = node_ref.borrow_mut();
    for i in 0..index {
        let mut basic = basis[i].borrow_mut();
        if basic.row == node.row || basic.col == node.col {
            basic.adjacency.push(Rc::clone(node_ref));
            node.adjacency.push(Rc::clone(&basis[i]));
        }
    }
}

fn remove_basic(basis: &mut Basis, index: usize, node_ref: &BasicRef) {
    // Find node in list
    let mut found = false;
    let node_id = node_ref.borrow().id;
    let mut node_index: usize = 0;
    for i in 0..index {
        if basis[i].borrow().id == node_id {
            found = true;
            node_index = i;
            break;
        }
    }
    assert!(found, "Node not found in basis.");

    for i in 0..node_ref.borrow().adjacency.len() {
        let adj_node_ref = Rc::clone(&node_ref.borrow().adjacency[i]);
        let mut adj_node = adj_node_ref.borrow_mut();
        
        let mut remove_index: Option<usize> = None;
        for j in 0..adj_node.adjacency.len() {
            if adj_node.adjacency[j].borrow().id == node_id {
                remove_index = Some(j);
            }
        }

        if let Some(remove_index_value) = remove_index {
            // Remove the node
            adj_node.adjacency.remove(remove_index_value);
            for i in remove_index_value..adj_node.adjacency.len() {
                adj_node.adjacency[i] = Rc::clone(&adj_node.adjacency[i + 1]);
            }
            adj_node.adjacency.pop();
        }
    }

    basis[node_index] = Rc::clone(&basis[index]);
}

fn reset_cursor(basis: &mut Basis, index: usize) {
    for i in 0..index {
        let mut basic = basis[i].borrow_mut();
        if basic.adjacency.len() > 0 {
            basic.cursor = Some(basic.adjacency.len() - 1);
        } else {
            basic.cursor = None;
        }
        basic.back = None;
        basic.color = Color::WHITE;
    }
}

pub fn emd(
    weight_x: &Vec<f64>,
    weight_y: &Vec<f64>,
    cost: &Vec<Vec<f64>>,
) -> f64 {
    let n_x = weight_x.len();
    let n_y = weight_y.len();

    let mut basis = initialize_flow(weight_x, weight_y, cost);

    // Iterate until optimality conditions satisfied
    let mut dual_x = vec![0.0; n_x];
    let mut dual_y = vec![0.0; n_y];
    let mut solved_x = vec![false; n_x];
    let mut solved_y = vec![false; n_y];

    let b = (n_x + n_y - 1) as usize;

    loop {
        for i in 0..n_x {
            solved_x[i] = false;
        }
        for i in 0..n_y {
            solved_y[i] = false;
        }
        reset_cursor(&mut basis, b);

        let mut var_ref = Rc::clone(&basis[0]);
        dual_x[var_ref.borrow().row] = 0.0;
        solved_x[var_ref.borrow().row] = true;

        loop {
            let update_var_ref: BasicRef;
            {
                let mut adj_cursor: Option<usize> = None;
                let mut var = var_ref.borrow_mut();

                var.color = Color::GRAY;
                if solved_x[var.row] {
                    dual_y[var.col] = cost[var.row][var.col] - dual_x[var.row];
                    solved_y[var.col] = true;
                } else if solved_y[var.col] {
                    dual_x[var.row] = cost[var.row][var.col] - dual_y[var.col];
                    solved_x[var.row] = true;
                } else {
                    panic!("Failed");
                }

                if let Some(cursor_value) = var.cursor {
                    for i in (0..=cursor_value).rev() {
                        let adj_var = var.adjacency[i].borrow();
                        if adj_var.color == Color::WHITE {
                            adj_cursor = Some(i);
                            break;
                        }
                    }
                }
                if let Some(adj_cursor_value) = adj_cursor {
                    if adj_cursor_value > 0 {
                        var.cursor = Some(adj_cursor_value - 1);
                    } else {
                        var.cursor = None;
                    }
                    let mut adj_var = var.adjacency[adj_cursor_value].borrow_mut();
                    adj_var.back = Some(Rc::clone(&var_ref));
                    update_var_ref = Rc::clone(&var.adjacency[adj_cursor_value]);
                } else {
                    var.color = Color::BLACK;

                    match var.back {
                        Some(ref back_value) => {
                            update_var_ref = Rc::clone(back_value);
                        },
                        None => {
                            break;
                        },
                    }
                }
            }
            var_ref = update_var_ref;
        }

        // Check for optimality
        let mut min_init = false;
        let mut min_row: usize = 0;
        let mut min_col: usize = 0;
        let mut min_slack = 0.0;

        for i in 0..n_x {
            for j in 0..n_y {
                let slack = cost[i][j] - dual_x[i] - dual_y[j];
                if !min_init || slack < min_slack {
                    min_init = true;
                    min_row = i;
                    min_col = j;
                    min_slack = slack;
                }
            }
        }

        let mut is_optimal = false;
        for i in 0..b {
            let basic = basis[i].borrow();
            // If the pivot variable is any of the
            // basis variables, then the optimal
            // solution has been found; set
            // min_slack = 0.0 explicitly to avoid
            if basic.row == min_row && basic.col == min_col {
                min_slack = 0.0;
                is_optimal = true;
                break;
            }
        }

        if min_slack >= -EPSILON || is_optimal {
            break;
        }

        // Introduce a new variable
        let mut new_var_ref = init_basic(min_row, min_col, 0.0);
        insert_basic(&mut basis, b, &new_var_ref);

        let root_ref: BasicRef = Rc::clone(&new_var_ref);
        let root_id = root_ref.borrow().id;
        reset_cursor(&mut basis, b + 1);

        loop {
            let update_var_ref: BasicRef;
            {
                let mut var = new_var_ref.borrow_mut();
                var.color = Color::GRAY;
                
                let mut adj_cursor: Option<usize> = None;
                if let Some(cursor_value) = var.cursor {
                    for i in (0..=cursor_value).rev() {
                        let adj_var = var.adjacency[i].borrow();
                        match adj_var.back {
                            Some(ref back_value) => {
                                let back_var = back_value.borrow();
                                if back_var.row == adj_var.row || back_var.col == adj_var.col {
                                    continue;
                                }
                            },
                            None => {},
                        }
                        if adj_var.id == root_id {
                            // Found a cycle
                            adj_cursor = Some(i);
                            break;
                        }
                        if adj_var.color == Color::WHITE {
                            adj_cursor = Some(i);
                            break;
                        }
                    }
                }
                if let Some(adj_cursor_value) = adj_cursor {
                    let is_gray: bool;
                    {
                        let adj_var = var.adjacency[adj_cursor_value].borrow();
                        is_gray = adj_var.color == Color::GRAY;
                    }
                    if is_gray {
                        // Found a cycle
                        let mut mut_root = root_ref.borrow_mut();
                        mut_root.back = Some(Rc::clone(&new_var_ref));
                        break;
                    } else {
                        if adj_cursor_value > 0 {
                            var.cursor = Some(adj_cursor_value - 1);
                        } else {
                            var.cursor = None;
                        }
                        let mut adj_var = var.adjacency[adj_cursor_value].borrow_mut();
                        adj_var.back = Some(Rc::clone(&new_var_ref));
                        update_var_ref = Rc::clone(&var.adjacency[adj_cursor_value]);
                    }
                } else {
                    var.color = Color::BLACK;

                    match var.back {
                        Some(ref back_value) => {
                            update_var_ref = Rc::clone(back_value);
                        },
                        None => {
                            panic!("Could not find a cycle");
                        },
                    }
                }
            }
            new_var_ref = update_var_ref;
        }

        // Find the largest flow that can be subtracted
        let mut sign = -1.0;
        let mut min_flow = 0.0;
        let mut to_remove: Option<BasicRef> = None;

        let mut loop_var_ref: BasicRef;
        match root_ref.borrow().back {
            Some(ref back_value) => {
                loop_var_ref = Rc::clone(back_value);
            },
            None => {
                panic!("Must have back value");
            },
        }

        loop {
            let update_var_ref: BasicRef;
            {
                let var = loop_var_ref.borrow();
                if sign < 0.0 && (to_remove.is_none() || var.flow < min_flow) {
                    min_flow = var.flow;
                    to_remove = Some(Rc::clone(&loop_var_ref));
                }
                sign *= -1.0;

                match var.back {
                    Some(ref back_value) => {
                        if back_value.borrow().id == root_id {
                            break;
                        }
                        update_var_ref = Rc::clone(back_value);
                    },
                    None => {
                        panic!("Must have back value");
                    },
                }
            }

            loop_var_ref = update_var_ref;
        }

        // Adjust flows
        sign = -1.0;
        root_ref.borrow_mut().flow = min_flow;

        match root_ref.borrow().back {
            Some(ref back_value) => {
                loop_var_ref = Rc::clone(back_value);
            },
            None => {
                panic!("Must have back value");
            },
        }

        loop {
            let update_var_ref: BasicRef;
            {
                let mut var = loop_var_ref.borrow_mut();
                var.flow += sign * min_flow;
                sign *= -1.0;

                match var.back {
                    Some(ref back_value) => {
                        if back_value.borrow().id == root_id {
                            break;
                        }
                        update_var_ref = Rc::clone(back_value);
                    },
                    None => {
                        panic!("Must have back value");
                    },
                }
            }
            
        }

        // Remove the basic variable that went to zero
        if let Some(var_to_remove) = to_remove {
            remove_basic(&mut basis, b, &var_to_remove);
        }
    }

    let mut distance = 0.0;
    for i in 0..b {
        let basic = basis[i].borrow();
        distance += basic.flow * cost[basic.row][basic.col];
    }

    distance
}

pub fn vectors_euclidean_distance(v1: &Vec<f64>, v2: &Vec<f64>) -> f64 {
    v1.iter()
      .zip(v2.iter())
      .map(|(x,y)| (x - y).powi(2))
      .sum::<f64>()
      .sqrt()
}

pub fn distance(v1: &Vec<f64>, v2: &Vec<f64>) -> f64 {
    let weight1 = vec![1./(v1.len() as f64); v1.len()];
    let weight2 = vec![1./(v2.len() as f64); v2.len()];

    let mut cost = Vec::with_capacity(v1.len());
    for i in v1.iter() {
        let mut cost_row = Vec::with_capacity(v2.len());
        for j in v2.iter() {
            cost_row.push(vectors_euclidean_distance(&vec![*i], &vec![*j]));
        }
        cost.push(cost_row)
    }

    emd(&weight1, &weight2, &cost)
}

#[cfg(test)]
mod tests {
    use crate::*;

    #[test]
    fn test_distance() {
        let x = vec![4., 3.];
        let y = vec![3., 4.];
        assert_eq!(distance(&x, &y), 0.);

        let x = vec![4., 3.];
        let y = vec![3., 5.];
        assert_eq!(distance(&x, &y), 0.5);

        let x = vec![4., 3.];
        let y = vec![3., 5., 3., 2.];
        assert_eq!(distance(&x, &y), 0.75);
    }
}
