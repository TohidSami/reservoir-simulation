# app.py
from flask import Flask, jsonify , request , send_file
import json
import time
from sqlalchemy.exc import OperationalError 
from run_model import map_db_to_engine , run_simulation
from models import db, SimulationCase, SimulationResult
import numpy as np
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import io
import os
def create_app():
    app = Flask(__name__)

    db_url= os.environ.get('DATABASE_URL')
    if not db_url:
       

        db_url='postgresql://admin:secretpassword@localhost:5433/reservoir_simulation'
    
    app.config['SQLALCHEMY_DATABASE_URI'] =db_url 
    app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False

    db.init_app(app)
    
    with app.app_context():
        retries=5
        while retries > 5:
           try:
            db.create_all()
            print("--- Database Tables Checked/Created Successfully ---")
            break
           except OperationalError:
              print(f"Database not ready yet... Retrying in 2 seconds ({retries} left)")
              time.sleep(2)
              retries -= 1
        if retries == 0:
            print("Error: Could not connect to Database after 5 attempts.")
    return app

app = create_app()

@app.route('/')
def home():
    return jsonify({
        "status": "Server is Running",
        "message": "Connected to PostgreSQL Database via Docker!"
    })


@app.route('/api/simulation' , methods=['POST'])
def run_simulation_api():
   try:
    data=request.json
    new_case = SimulationCase(
            name=data.get('name', 'Unnamed Case'),
            status='running',
            nx=data.get('nx', 15),
            ny=data.get('ny', 15),
            nz=data.get('nz', 3),
            swc=data.get('swc', 0.1),
            sor=data.get('sor', 0.3),
            inj_rate_target=data.get('inj_rate_target', 1500.0),
            prod_bhp_min=data.get('prod_bhp_min', 3000.0),
            dt=data.get('dt', 10.0),
            n_steps=data.get('n_steps', 20) )
    db.session.add(new_case)
    db.session.commit()
    print(f"Case saved with ID: {new_case.id}")
    case_engine=map_db_to_engine(new_case)
    sw_t,p_t,time_t,well=run_simulation(case_engine)
    for i,time_val in enumerate(time_t):
        sw_json=json.dumps(sw_t[i].tolist())
        p_json=json.dumps(p_t[i].tolist())
        current_well={
           'q_inj':well['q_inj_w'][i],
           "q_oil": well['q_prod_o'][i],
            "q_wat": well['q_prod_w'][i],
            "bhp_inj": well['pbh_inj'][i],
            "bhp_prod": well['pbh_prod'][i]

        }
        well_json=json.dumps(current_well)

        result=SimulationResult(
        case_id=new_case.id,
        step_number = i+1,
        time_elapsed=time_val,
        saturation_data=sw_json,
        pressure_data=p_json,
        well_dynamic_data=well_json
        )
        db.session.add(result)
    new_case.status='completed'
    db.session.commit()
    return jsonify({"message":"simulation completed successfully",
                    "case_id":new_case.id,
                    "step":len(time_t),
                    "status":"completed"}),201
   
   except Exception as e:
      return jsonify ({"error ":str(e)}),500


@app.route("/api/simulation/<int:case_id>",methods=["GET"])
def get_results(case_id):
   try:
    case=SimulationCase.query.get(case_id)
    if not case:
        return jsonify({"error": "Case not found"}), 404
   

    result=SimulationResult.query.filter_by(case_id=case_id).order_by(SimulationResult.step_number).all()
    output={
        "case_id":case.id,
        "case_name":case.name,
        "status":case.status,
        "grid":{"nx":case.nx , "ny":case.ny,"nz":case.nz},
        "steps":[]
    }
    for res in result:
        p_mat=json.loads(res.pressure_data) 
        sw_mat=json.loads(res.saturation_data)
        step_info={
            "step":res.step_number,
            "time_days":res.time_elapsed,
            "pressure_map":p_mat,
            "saturation_map":sw_mat,
        }
        output["steps"].append(step_info)
    return jsonify (output),200
   except Exception as e:
        print(f"Error fetching results: {e}")
        return jsonify({"error": str(e)}), 500
   
@app.route("/api/simulation/<int:case_id>/plot",methods=['GET'])
def plot_meanP(case_id):
   try:
    result=SimulationResult.query.filter_by(case_id=case_id).order_by(SimulationResult.step_number).all()
    if not result: return jsonify({"error":" not found data "}),404
    
    p=[]
    t=[]
    for res in result:
       p_mat=json.loads(res.pressure_data)
       t_mat=res.time_elapsed
       p.append(np.nanmean(p_mat))
       t.append(t_mat)
    fig,ax=plt.subplots(figsize=(10, 6))   
    # plt.figure(figsize=(10,6))
    plt.plot(t,p, marker='o', linestyle='-', color='b', label='Field Avg Pressure')
    plt.title(f"Simulation Result: {case_id}")
    plt.xlabel("Time (Days)")
    plt.ylabel("Average Reservoir Pressure (psi)")
    plt.grid(True)
    plt.legend()
    
    print("Displaying Plot...")
    # plt.show()
    buf=io.BytesIO()
    plt.savefig(buf,format='png')
    buf.seek(0)
    plt.close(fig)
    return send_file(buf, mimetype='image/png')
   except Exception as e:
      print(f"Plot Error: {e}")  
      return jsonify({"error": str(e)}), 500


@app.route("/api/simulation/<int:case_id>/well_plot")
def plotWell(case_id):
     try: 
        result=SimulationResult.query.filter_by(case_id=case_id).order_by(SimulationResult.step_number).all()
        if not result: return jsonify({"error":"404 not found this case"}),404
        t=[]
        q_oil=[]
        q_water=[]
        q_inj=[]
        bhp_p=[]
        bhp_i=[]
        for res in result:
           t.append(res.time_elapsed)
           if res.well_dynamic_data:
                w_data=json.loads(res.well_dynamic_data)
                q_oil.append(w_data.get('q_oil',0))
                q_inj.append(w_data.get('q_inj',0))
                q_water.append(w_data.get('q_wat',0))
                bhp_i.append(w_data.get('bhp_inj',0))
                bhp_p.append(w_data.get('bhp_prod',0))
           else:
              q_oil.append(0); q_water.append(0); q_inj.append(0)
              bhp_p.append(0); bhp_i.append(0)
            
        fig,ax1=plt.subplots(figsize=(12,7))
        l1=ax1.plot(t,q_oil,color='green',linestyle='-',linewidth=2,label='Oil Rate')
        l2=ax1.plot(t,q_water,color='blue',linestyle='--',linewidth=2,label='Water Rate')
        l3=ax1.plot(t,q_inj,color='cyan',linestyle='-.',linewidth=2,label='Water injection Rate')
        ax1.set_xlabel('Time (dayes)',fontsize=12)
        ax1.set_ylabel('Rate (STB/day)',fontsize=12)
        ax1.grid(True,linestyle=':' , alpha=0.5)

        ax2=ax1.twinx()
        l4=ax2.plot(t,bhp_p,color='red',linestyle='-',marker='o',markersize=4, label='BHP Producer')
        l5 = ax2.plot(t, bhp_i, color='purple', linestyle='-', marker='x', 
        markersize=4, label='BHP Injector')
        ax2.set_ylabel('Pressure (psi)', fontsize=12)
        lines = l1 + l2 + l3 + l4 + l5
        labels = [l.get_label() for l in lines]
        ax1.legend(lines, labels, loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=5, frameon=True)
        
        plt.title(f"Production & Injection Profile - Case {case_id}", fontsize=14, pad=40)
        plt.tight_layout()

        
        buf = io.BytesIO()
        plt.savefig(buf, format='png', dpi=100)
        buf.seek(0)
        plt.close(fig)
        
        return send_file(buf, mimetype='image/png')
     except Exception as e:
        return jsonify({f"plot error": str(e)})   





if __name__ == '__main__':
    app.run(debug=True, port=5000 , host='0.0.0.0')