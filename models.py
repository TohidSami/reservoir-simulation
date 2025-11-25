# models.py
from flask_sqlalchemy import SQLAlchemy
from datetime import datetime


db = SQLAlchemy()


class SimulationCase(db.Model):
    __tablename__ = 'simulation_cases'  

    id = db.Column(db.Integer, primary_key=True)   
    name = db.Column(db.String(100), nullable=False) 
    created_at = db.Column(db.DateTime, default=datetime.utcnow) 
    status = db.Column(db.String(20), default='pending') # 

    nx = db.Column(db.Integer, default=15)
    ny = db.Column(db.Integer, default=15)
    nz = db.Column(db.Integer, default=3)

    swc = db.Column(db.Float, default=0.1) 
    sor = db.Column(db.Float, default=0.3)  
    
    pc_coeff = db.Column(db.Float, default=5.0)
    pc_exp = db.Column(db.Float, default=2.0)
    krw_max = db.Column(db.Float, default=0.4)
    krw_exp = db.Column(db.Float, default=1.2)
    kro_max = db.Column(db.Float, default=0.35)
    kro_exp = db.Column(db.Float, default=2.5)

    bo_slope = db.Column(db.Float, default=-6e-6)
    bo_intercept = db.Column(db.Float, default=1.275)
    mu_oil_low = db.Column(db.Float, default=0.5)
    mu_oil_high = db.Column(db.Float, default=0.6)
    p_low = db.Column(db.Float, default=2500.0)
    p_high = db.Column(db.Float, default=3500.0)

    bw_base = db.Column(db.Float, default=1.03)
    cw = db.Column(db.Float, default=3e-6)
    muw = db.Column(db.Float, default=0.4)
    p_ref = db.Column(db.Float, default=4000.0)

    inj_rate_target = db.Column(db.Float, default=1500.0)
    inj_bhp_max = db.Column(db.Float, default=6000.0)
    prod_rate_target = db.Column(db.Float, default=300.0)
    prod_bhp_min = db.Column(db.Float, default=3000.0)
    dt = db.Column(db.Float, default=10.0)
    n_steps = db.Column(db.Integer, default=50)
    scenario_mode = db.Column(db.Integer, default=3)

    results = db.relationship('SimulationResult', backref='case', lazy=True)

    def __repr__(self):
        return f"<Case {self.id}: {self.name}>"


class SimulationResult(db.Model):
    __tablename__ = 'simulation_results'

    id = db.Column(db.Integer, primary_key=True)

    case_id = db.Column(db.Integer, db.ForeignKey('simulation_cases.id'), nullable=False)

    step_number = db.Column(db.Integer, nullable=False)
    time_elapsed = db.Column(db.Float)

    pressure_data = db.Column(db.Text)      
    saturation_data = db.Column(db.Text)    

    well_dynamic_data = db.Column(db.Text)

    def __repr__(self):
        return f"<Result Case:{self.case_id} Step:{self.step_number}>"