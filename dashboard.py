import streamlit as st
import requests
from PIL import Image
import io

# تنظیمات صفحه
st.set_page_config(page_title="Reservoir Sim Dashboard", page_icon="🛢️", layout="wide")

# آدرس API (چون داشبورد روی ویندوز اجرا می‌شود، باید به لوکال‌هاست وصل شود)
BASE_URL = "http://localhost:5000/api/simulation"

st.title("🛢️ Reservoir Simulation Dashboard")
st.markdown("---")

# ==========================================
# 1. بخش ورودی‌ها (Sidebar)
# ==========================================
with st.sidebar:
    st.header("⚙️ Simulation Parameters")
    
    # نام کیس
    case_name = st.text_input("Case Name", value="Test_Case_01")
    
    # گرید
    st.subheader("Grid Dimensions")
    c1, c2, c3 = st.columns(3)
    nx = c1.number_input("NX", min_value=5, max_value=100, value=15)
    ny = c2.number_input("NY", min_value=5, max_value=100, value=15)
    nz = c3.number_input("NZ", min_value=1, max_value=20, value=3)
    
    # خواص سنگ و سیال
    st.subheader("Fluid Properties")
    swc = st.slider("Connate Water (Swc)", 0.0, 0.5, 0.1)
    sor = st.slider("Residual Oil (Sor)", 0.0, 0.5, 0.3)
    
    # چاه‌ها
    st.subheader("Well Controls")
    inj_rate = st.number_input("Injection Rate (STB/D)", value=1500.0)
    prod_rate = st.number_input("Production Rate (STB/D)", value=500.0)
    prod_bhp = st.number_input("Min Prod BHP (psi)", value=3000.0)
    
    # زمان
    st.subheader("Time Settings")
    dt = st.number_input("Delta T (Days)", value=10.0)
    n_steps = st.number_input("Number of Steps", value=20)

# ==========================================
# 2. بخش اصلی (Tabs)
# ==========================================
tab1, tab2 = st.tabs(["🚀 Run Simulation", "📊 View Plots"])

# --- Tab 1: اجرای شبیه سازی ---
with tab1:
    st.info("👈 Set your parameters in the sidebar, then click 'Run Simulation'.")
    
    if st.button("Run Simulation", type="primary"):
        # آماده‌سازی داده برای ارسال (JSON Payload)
        payload = {
            "name": case_name,
            "nx": nx, "ny": ny, "nz": nz,
            "swc": swc, "sor": sor,
            "inj_rate_target": inj_rate,
            "prod_rate_target": prod_rate, # اضافه شده طبق خواسته شما
            "prod_bhp_min": prod_bhp,
            "dt": dt,
            "n_steps": n_steps
        }
        
        # نمایش وضعیت در حال اجرا
        with st.spinner("Connecting to Docker API and running simulation..."):
            try:
                # درخواست POST به API
                response = requests.post(BASE_URL, json=payload)
                
                if response.status_code == 201:
                    data = response.json()
                    st.success(f"✅ Simulation Completed Successfully!")
                    
                    # نمایش اطلاعات ذخیره شده
                    st.json(data)
                    
                    # ذخیره ID برای استفاده راحت‌تر در تب بعدی
                    st.session_state['last_case_id'] = data['case_id']
                else:
                    st.error(f"❌ Error: {response.text}")
                    
            except requests.exceptions.ConnectionError:
                st.error("❌ Connection Failed! Is Docker running? (http://localhost:5000)")

# --- Tab 2: نمایش پلات‌ها ---
with tab2:
    st.write("Enter a Case ID to view its results.")
    
    # اگر قبلاً کیسی ران شده، آیدی آن را پیش‌فرض بگذار
    default_id = st.session_state.get('last_case_id', 1)
    case_id_input = st.number_input("Case ID", min_value=1, value=default_id, step=1)
    
    col_plot1, col_plot2 = st.columns(2)
    
    # دکمه اول: پلات فشار میانگین
    with col_plot1:
        if st.button("📈 Show Average Pressure Plot"):
            with st.spinner("Fetching plot..."):
                try:
                    # درخواست GET به API پلات فشار
                    plot_url = f"{BASE_URL}/{case_id_input}/plot"
                    resp = requests.get(plot_url)
                    
                    if resp.status_code == 200:
                        # تبدیل بایت‌های تصویر به عکس قابل نمایش
                        image = Image.open(io.BytesIO(resp.content))
                        st.image(image, caption=f"Average Pressure - Case {case_id_input}", use_column_width=True)
                    else:
                        st.error("Plot not found. Check Case ID.")
                except Exception as e:
                    st.error(f"Error: {e}")

    # دکمه دوم: پلات چاه
    with col_plot2:
        if st.button("📉 Show Well Production Plot"):
            with st.spinner("Fetching plot..."):
                try:
                    # درخواست GET به API پلات چاه
                    well_url = f"{BASE_URL}/{case_id_input}/well_plot"
                    resp = requests.get(well_url)
                    
                    if resp.status_code == 200:
                        image = Image.open(io.BytesIO(resp.content))
                        st.image(image, caption=f"Well Profile - Case {case_id_input}", use_column_width=True)
                    else:
                        st.error("Well data not found.")
                except Exception as e:
                    st.error(f"Error: {e}")